"""Publish a new immutable-version release using the workflow's GITHUB_TOKEN.

Only a new version or an incomplete draft at the same commit can be published.
Existing tags/releases are never retargeted or deleted.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import quote
from urllib.request import Request, urlopen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', default='release_dist')
    a = parser.parse_args()
    root = Path(a.directory)
    meta = json.loads((root / 'release_provenance.json').read_text())
    sha, version = meta['commit'], meta['version']
    repository = os.environ['GITHUB_REPOSITORY']
    if repository != 'Shiywa/conserved-host-response-framework':
        raise ValueError('Publication is restricted to this repository')
    if sha != os.environ['GITHUB_SHA']:
        raise ValueError('Archive commit differs from the verified workflow commit')
    token = os.environ['GITHUB_TOKEN']
    base = 'https://api.github.com/repos/' + repository

    def api(path, payload=None, *, method='GET', missing=False, raw=False):
        data = payload if raw else (json.dumps(payload).encode() if payload is not None else None)
        headers = {'Authorization': 'Bearer ' + token, 'User-Agent': 'CRRM388-release',
                   'Accept': 'application/vnd.github+json', 'X-GitHub-Api-Version': '2022-11-28',
                   'Content-Type': 'application/octet-stream' if raw else 'application/json'}
        url = path if path.startswith('https://uploads.github.com/') else base + path
        try:
            with urlopen(Request(url, data=data, headers=headers, method=method), timeout=120) as response:
                return json.load(response)
        except HTTPError as error:
            if missing and error.code == 404:
                return None
            # Never serialize authorization headers or environment variables.
            raise RuntimeError(f'GitHub {method} failed with HTTP {error.code}') from None

    tag = 'v' + version
    release = api('/releases/tags/' + tag, missing=True)
    ref = api('/git/ref/tags/' + tag, missing=True)
    if ref:
        obj = ref['object']
        while obj['type'] == 'tag':
            obj = api('/git/tags/' + obj['sha'])['object']
        if obj['sha'] != sha:
            raise ValueError('Version already points to another commit; choose a new version')
    if release and not release['draft']:
        raise ValueError('Published release is protected; it will not be modified')
    if not ref:
        api('/git/refs', {'ref': 'refs/tags/' + tag, 'sha': sha}, method='POST')
    body = (f'Current CRRM388 manuscript workflows (388 genes; 280 up / 108 down).\n\n'
            f'Commit: `{sha}`\n\n'
            'Code-only release: signed fit-reference/project scoring, Results 1–6 analysis and figure code, '
            'Tables S1–S8 generation/export, and support for local-input clinical replay. Matrices and complete result tables are not distributed. '
            'See LINUX_VALIDATION_REPORT.md for executed tests and scope. '
            'Raw DE and the full large-cell resource were not freshly rerun. '
            'The v0.1.0 release/tag remain unchanged.\n\n'
            'Verify the source archive against SHA256SUMS; the archive contains RELEASE_PROVENANCE.json.')
    if not release:
        release = api('/releases', dict(tag_name=tag, target_commitish=sha,
                                       name=tag + ' — current CRRM388 study', body=body,
                                       draft=True, prerelease=False), method='POST')
    existing = {x['name']: x for x in release.get('assets', [])}
    names = [f'conserved-host-response-framework-{version}.tar.gz', 'SHA256SUMS']
    for name in names:
        content = (root / name).read_bytes()
        if name in existing:
            digest = existing[name].get('digest')
            if digest != 'sha256:' + hashlib.sha256(content).hexdigest():
                raise ValueError('Existing draft asset cannot be verified; no overwrite attempted')
            continue
        api(release['upload_url'].split('{')[0] + '?name=' + quote(name), content, method='POST', raw=True)
    published = api('/releases/' + str(release['id']), dict(draft=False), method='PATCH')
    print(json.dumps(dict(tag=tag, commit=sha, release_url=published['html_url'], published=True)))


if __name__ == '__main__':
    main()
