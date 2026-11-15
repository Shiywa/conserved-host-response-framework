"""Create a deterministic archive from one committed Git tree, never the cwd."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
import tomllib


def git(*arguments):
    return subprocess.check_output(['git', *arguments])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ref', default='HEAD')
    p.add_argument('--output-dir', default='release_dist')
    a = p.parse_args()
    sha = git('rev-parse', '--verify', a.ref + '^{commit}').decode().strip()
    tree = git('rev-parse', sha + '^{tree}').decode().strip()
    version = tomllib.loads(git('show', sha + ':pyproject.toml').decode())['project']['version']
    prefix = f'conserved-host-response-framework-{version}'
    raw = git('archive', '--format=tar', '--prefix=' + prefix + '/', sha)
    timestamp = int(git('show', '-s', '--format=%ct', sha))
    buffer = io.BytesIO(raw)
    provenance = json.dumps(dict(version=version, commit=sha, tree=tree,
                                repository='https://github.com/Shiywa/conserved-host-response-framework'), indent=2).encode() + b'\n'
    with tarfile.open(fileobj=buffer, mode='a') as tar:
        entry = tarfile.TarInfo(prefix + '/RELEASE_PROVENANCE.json')
        entry.size = len(provenance); entry.mtime = timestamp; entry.mode = 0o644
        tar.addfile(entry, io.BytesIO(provenance))
    out = Path(a.output_dir); out.mkdir(parents=True, exist_ok=True)
    archive = out / (prefix + '.tar.gz')
    with archive.open('wb') as f, gzip.GzipFile(filename='', mode='wb', fileobj=f, mtime=0) as z:
        z.write(buffer.getvalue())
    checksum = hashlib.sha256(archive.read_bytes()).hexdigest()
    (out / 'SHA256SUMS').write_text(checksum + '  ' + archive.name + '\n')
    (out / 'release_provenance.json').write_bytes(provenance)
    print(json.dumps(dict(commit=sha, version=version, archive=archive.name, sha256=checksum)))


if __name__ == '__main__':
    main()
