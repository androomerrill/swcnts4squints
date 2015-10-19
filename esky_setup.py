from __future__ import unicode_literals
import sys

from esky import bdist_esky
from esky.bdist_esky import Executable
import ez_setup
ez_setup.use_setuptools()
#from distutils.core import setup
from setuptools import setup, find_packages

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
            
    extra_options = dict(setup_requires=['py2exe'],
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
                         data_files=[("phonon_backend", 
                        ["C:\Python27\Lib\site-packages\PyQt4\plugins\phonon_backend\phonon_ds94.dll"])])

    setup(name='swcnts4squints',
          version=0.5,
          description="python scripts",
          long_description="""\
                  Simple app for calculating physical properties of carbon nanotubes""",
          classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
          keywords='python scripts',
          author='Andrew Merril',
          author_email='androomerrill@gmail.com',
          url='https://github.com/androomerrill/swcnts4squints',
          license='BSD 2-Clause',
          packages=find_packages(exclude=['swcnts4squints']),
          #package_dir={'': 'swcnts4squints'},
          include_package_data=True,
          #exclude_package_data={'': ['README.md']},
          #zip_safe=False,
          install_requires=['numpy', 'scipy'],
          entry_points={
              'gui_scripts': [
                  'swcnts4squints = swcnts4squints:main',
                  ],
            },
          **extra_options
    )
    
    for f in py_files:
        os.remove(f)
                         
else:
    extra_options = {}
    if sys.platform == 'darwin':
        #extra_options = dict(setup_requires=['py2app'],
        extra_options = dict(app=['swcnts4squints/swcnts4squints.py'],
                             scripts=[Executable('swcnts4squints/swcnts4squints.py',)],
                             options={'bdist_esky': {'freezer_module': 'py2app',
                                                     'freezer_options': 
                                                        {"argv_emulation": True,
                                                         "iconfile": 'images/swcnts4squints.icns',
                                                         "includes": 
                                                            ['sip', 'PyQt4.QtCore', 'PyQt4.QtGui']}}})
                                            
    
    setup(name='swcnts4squints',
          version=0.5,
          description="python scripts",
          long_description="""\
                  Simple app for calculating physical properties of carbon nanotubes""",
          classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
          keywords='python scripts',
          author='Andrew Merril',
          author_email='androomerrill@gmail.com',
          url='https://github.com/androomerrill/swcnts4squints',
          license='BSD 2-Clause',
          packages=find_packages(),
          #package_dir={'': 'swcnts4squints'},
          include_package_data=True,
          #exclude_package_data={'': ['README.md']},
          zip_safe=False,
          install_requires=['numpy', 'scipy'],
          #entry_points={
          #    'gui_scripts': [
          #        'swcnts4squints = swcnts4squints.swcnts4squints:main',
          #        ],
          #},
          **extra_options
    )

