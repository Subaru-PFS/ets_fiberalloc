from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import subprocess

__version__ = "2.1"

try:
    tmp = subprocess.run(["git", "describe", "--dirty"], capture_output=True, text=True)
    gitversion = tmp.stdout
except:
    gitversion = "unknown.gitversion"

__version__ = __version__ + "+" + gitversion

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

extsrc = ['src/external/string_utils.cc', 'src/external/error_handling.cc']

ext_modules = [
    Extension(
        'pyETS',
        ['src/pyETS.cc','src/ets.cc'] + extsrc,
        include_dirs=['src','src/external',
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True)
        ],
        language='c++'
    ),
    ]

# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.
    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)

setup(name="ets_fiber_assigner",
      version=__version__,
      author="Martin Reinecke",
      author_email="martin@mpa-garching.mpg.de",
      description="Fiber assignment code for the PFS instrument",
      packages=find_packages(include=["ets_fiber_assigner", "ets_fiber_assigner.*"]),
      zip_safe=False,
      cmdclass={'build_ext': BuildExt},
      ext_modules=ext_modules,
      dependency_links = [],
      extras_require={'pybind11': 'pybind11>=2.2.1'},
      install_requires = ["pybind11>=2.2.1"],
      license="GPLv2",
)

