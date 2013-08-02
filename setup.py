import sys
major, minor1, minor2, s, tmp = sys.version_info

#if major==2 and minor1 < 7 or major < 2:
#    raise SystemExit("""cntapp requires Python 2.6 or later""")

if sys.platform == 'darwin':
    from setuptools import setup

    APP = ['swcnts4squints.py']
    OPTIONS = {'argv_emulation': True,
               'includes': ['sip', 'PyQt4.QtCore', 'PyQt4.QtGui']}

    setup(
        app=APP,
        options={'py2app': OPTIONS},
        setup_requires=['py2app'],
    )

elif sys.platform == 'win32':
    from distutils.core import setup
    import py2exe

    setup(console=['swcnts4squints.py'],
          options={
              'py2exe': {
                  "bundle_files": 1,
                  "dll_excludes": [
                      "MSVCP90.dll",
                      "MSWSOCK.dll",
                      "mswsock.dll",
                      "powrprof.dll"],
                  "includes": [
                      'sip',
                      'PyQt4.QtCore',
                      'PyQt4.QtGui'],
              },
          },
          data_files=[('phonon_backend',
              ['C:\Python27\Lib\site-packages\PyQt4\plugins\phonon_backend\phonon_ds94.dll'])]
          )
