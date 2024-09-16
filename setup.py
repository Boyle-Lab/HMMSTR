from setuptools import setup, Extension
#from Cython.Build import cythonize
import subprocess
import numpy as np

#Custom function to find install location of glib-2.0
def pkgconfig(*packages):
    cmd = ['pkg-config', '--cflags', '--libs'] + list(packages)
    output = subprocess.check_output(cmd).decode().split()
    cflags = [flag[2:] for flag in output if flag.startswith('-I')]
    libs = [flag[2:] for flag in output if flag.startswith('-l')]
    return cflags, libs

#added [glib_cflags[0]+"/glib"] to designate where the included/hiddedn header files are for glib-2.0, TODO make sure this works for the install I include in README
glib_cflags, glib_libs = pkgconfig('glib-2.0')

numpy_include_dir = np.get_include()

try:
    from Cython.Build import cythonize
    extensions = cythonize([
    Extension(
        "c_files.hmm",
        ["src/c_files/hmm.pyx","src/c_files/sequence.c","src/c_files/viterbi.c","src/c_files/hmmutils.c","src/c_files/nrutil.c","src/c_files/hmmrand.c"],
        libraries=glib_libs,
        include_dirs=["src/c_files"]+glib_cflags+[glib_cflags[0]+"/glib"],
        extra_compile_args=[],
        extra_link_args=[],
    ),
    Extension(
            "pylibsais",
            sources=["src/c_files/libsais.c", "src/c_files/pylibsais.c"],
            include_dirs=[numpy_include_dir, numpy_include_dir + "/numpy", "src/c_files"],
            extra_compile_args=[],
            extra_link_args=[],
    )
])
except ImportError:
    extensions = [
    Extension(
        "c_files.hmm",
        ["src/c_files/hmm.c","src/c_files/sequence.c","src/c_files/viterbi.c","src/c_files/hmmutils.c","src/c_files/nrutil.c","src/c_files/hmmrand.c"],
        libraries=glib_libs,
        include_dirs=["src/c_files"]+glib_cflags+[glib_cflags[0]+"/glib"],
        extra_compile_args=[],
        extra_link_args=[],
    ),
    Extension(
        "pylibsais",
        sources=["src/c_files/libsais.c", "src/c_files/pylibsais.c"],
        include_dirs=[numpy_include_dir, numpy_include_dir + "/numpy", "src/c_files"],
        extra_compile_args=[],
        extra_link_args=[],
    )
]


setup(
    name="HMMSTR",
    version="1.0.5",
    python_requires='>=3.8.17',
    packages=["HMMSTR","GMM_stats","HMMSTR_utils","c_files","process_vit_res","profile_HMM","process_read","KDE_stats","motif_comp"],
    include_dirs=[np.get_include()],
    ext_modules=extensions,
    install_requires=[
        'colorama','numpy','pandas','pickleshare','scikit-learn','scipy','seaborn','importlib-resources','mappy','pysam','biopython',
        'matplotlib','multiprocess','umap-learn','levenshtein','pyabpoa'
    ],
    zip_safe=False,
    scripts=['src/c_files/test_hmm_cython.py'],
    # package_data = {
    #     'c_files': ['*.pxd']},
)