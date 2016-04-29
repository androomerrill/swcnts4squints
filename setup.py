"""SWCNTS4SQUINTS: interactive app for calculating nanotube properties.

"""
from __future__ import unicode_literals
from builtins import str

DOCLINES = __doc__.split("\n")

from setuptools import setup, find_packages
import os
import sys
import subprocess
import ez_setup
ez_setup.use_setuptools()

if sys.version_info[:2] < (2, 7):
    raise RuntimeError("Python version 2.7 required.")

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
License :: OSI Approved :: BSD License
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Programming Language :: Python
Programming Language :: Python :: 2.7
Topic :: Software Development
Topic :: Scientific/Engineering

"""

MAJOR = 0
MINOR = 5
MICRO = 1
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


def git_version():
    """Return the GIT version as a string."""
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = \
            subprocess.Popen(cmd,
                             stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


def get_version_info():
    # Adding the git rev number needs to be done inside
    # write_version_py(), otherwise the import of swcnts4squints.version messes
    # up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('swcnts4squints/version.py'):
        # must be a source distribution, use existing version file
        # load it as a separate module to not load swcnts4squints/__init__.py
        import imp
        version = imp.load_source('swcnts4squints.version',
                                  'swcnts4squints/version.py')
        GIT_REVISION = version.git_revision
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev-' + GIT_REVISION[:7]

    return FULLVERSION, GIT_REVISION


def write_version_py(filename='swcnts4squints/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM SWCNTS4SQUINTS SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    FULLVERSION, GIT_REVISION = get_version_info()

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()


def setup_package():
    # Rewrite the version file everytime
    write_version_py()

    FULLVERSION, GIT_REVISION = get_version_info()

    setup_options = \
        dict(name='swcnts4squints',
             version=FULLVERSION,
             author='Andrew Merrill',
             author_email='androomerrill@gmail.com',
             description=DOCLINES[0],
             long_description="\n".join(DOCLINES[2:]),
             url='https://github.com/androomerrill/swcnts4squints',
             license='BSD 2-Clause',
             classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
             platforms=["Windows", "Linux", "OS-X", "Unix"],
             test_suite='nose.collector',
             packages=find_packages(exclude=['docs', 'ez_setup', 'examples',
                                             'tests']),
             #package_data={'swcnts4squints': []},
             #package_dir = {'': 'swcnts4squints'},
             include_package_data=True,
             exclude_package_data={'':
                 ['README', 'README.rst', '*.gif', '*.html', '*.ui']},
             zip_safe=False,
             install_requires=['numpy>=1.7', 'scipy>=0.12',
                               'scikit-nano==0.2.25'],
             )

    extra_options = {}
    if sys.platform == 'win32':
        import os
        import py2exe
        import shutil
        pkg_dir = 'swcnts4squints'
        pkg_files = os.listdir(pkg_dir)
        py_files = []
        for f in pkg_files:
            fpath = os.path.join('swcnts4squints', f)
            if f.endswith('.py'):
                shutil.copyfile(fpath, f)
                py_files.append(f)

        extra_options = \
            dict(setup_requires=['py2exe'],
                 console=['swcnts4squints.py'],
                 options={'py2exe': {
                            "bundle_files": 1,
                            "dll_excludes": ["MSVCP90.dll",
                                             "MSWSOCK.dll",
                                             "mswsock.dll",
                                             "powrprof.dll"],
                            "includes": ["sip", "PyQt4.QtCore",
                                         "PyQt4.QtGui"]}},
                 windows=[{"script": "swcnts4squints.py",
                           "icon_resources":
                                [(1, "images/swcnts4squints.ico")]}],
                 data_files=[
                     ("phonon_backend", [
                         "C:\Python27\Lib\site-packages\PyQt4\plugins"
                         "\phonon_backend\phonon_ds94.dll"])],
                 entry_points={'gui_scripts':
                               ['swcnts4squints = swcnts4squints:main']},
                 )

        for f in py_files:
            os.remove(f)

    else:
        if sys.platform == 'darwin':
            extra_options = \
                dict(setup_requires=['py2app'],
                     app=['swcnts4squints/swcnts4squints.py'],
                     options={'py2app': {"argv_emulation": True,
                                         "iconfile":
                                         'images/swcnts4squints.icns',
                                         "includes":
                                         ['sip', 'PyQt4.QtCore',
                                          'PyQt4.QtGui']}})

        extra_options.update(
            dict(entry_points={
                 'gui_scripts':
                 ['swcnts4squints = swcnts4squints.swcnts4squints:main']}))

    setup_options.update(extra_options)
    setup(**setup_options)

if __name__ == '__main__':
    setup_package()
