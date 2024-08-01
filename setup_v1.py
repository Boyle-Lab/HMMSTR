from setuptools import setup, Extension, find_packages
from setuptools.command.install import install
import subprocess
import os

class CustomInstallCommand(install):
    def run(self):
        install_motifs()
        install.run(self)


def pkgconfig(*packages):
    cmd = ['pkg-config', '--cflags', '--libs'] + list(packages)
    output = subprocess.check_output(cmd).decode().split()
    cflags = [flag[2:] for flag in output if flag.startswith('-I')]
    libs = [flag[2:] for flag in output if flag.startswith('-l')]
    return cflags, libs

glib_cflags, glib_libs = pkgconfig('glib-2.0')

try:
    from Cython.Build import cythonize
    extensions = cythonize([
        Extension(
            "c_files.hmm",
            ["src/c_files/hmm.pyx", "src/c_files/sequence.c", "src/c_files/viterbi.c",
             "src/c_files/hmmutils.c", "src/c_files/nrutil.c", "src/c_files/hmmrand.c"],
            libraries=glib_libs,
            include_dirs=["src/c_files"] + glib_cflags + [glib_cflags[0] + "/glib"],
            extra_compile_args=[],
            extra_link_args=[],
        ),
    ])
except ImportError:
    extensions = [
        Extension(
            "c_files.hmm",
            ["src/c_files/hmm.c", "src/c_files/sequence.c", "src/c_files/viterbi.c",
             "src/c_files/hmmutils.c", "src/c_files/nrutil.c", "src/c_files/hmmrand.c"],
            libraries=glib_libs,
            include_dirs=["src/c_files"] + glib_cflags + [glib_cflags[0] + "/glib"],
            extra_compile_args=[],
            extra_link_args=[],
        ),
    ]

setup(
    name="HMMSTR",
    version="1.0.0",
    python_requires='>=3.8.17',
    packages=find_packages(),
    ext_modules=extensions,
    install_requires=['colorama','numpy','pandas','pickleshare','scikit-learn','scipy','seaborn','importlib-resources','mappy','biopython'],
    zip_safe=False,
    scripts=['src/c_files/test_hmm_cython.py'],
    cmdclass={
        'install': CustomInstallCommand,
    },
    # package_data = {
    #     'c_files': ['*.pxd']},
)

def install_motifs():
    conda_path = os.environ.get('CONDA_EXE', None)
    if not conda_path:
        raise EnvironmentError("Conda executable not found. Please ensure conda is installed and available in the PATH.")
    conda_base = os.path.dirname(os.path.dirname(conda_path))

    repo_url = "https://github.com/holstegelab/MotifScope.git"
    clone_cmd = ["git", "clone", repo_url]
    subprocess.check_call(clone_cmd)

    os.chdir("MotifScope")
    env_name = "motifscope"

    pylibsais_install_cmd = f"""
    source {conda_base}/etc/profile.d/conda.sh && \
    conda activate {env_name} && \
    sh INSTALL.sh
    """
    subprocess.run(pylibsais_install_cmd, shell=True, executable='/bin/bash') 

    motifscope_setup_cmd = f"""
    source {conda_base}/etc/profile.d/conda.sh && \
    conda activate {env_name} && \
    python setup.py install
    """
    subprocess.run(motifscope_setup_cmd, shell=True, executable='/bin/bash') 