
from setuptools import setup, find_packages
setup(
  name = 'boson_ladder',         # How you named your package folder (MyLib)
  packages = find_packages(),   # Chose the same as "name"
  version = '0.0.1',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Bosonic ladder operator calculator.',   # Give a short description about your library
  author = 'Hendry Minfui Lim',                   # Type in your name
  author_email = 'hendryadi01@gmail.com',      # Type in your E-Mail
  url = 'https://github.com/hendry24/boson_ladder',   # Provide either the link to your github or to your website
  download_url = '',    # I explain this later on
  keywords = ['Bosonic ladder operators'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'sympy>=1.13.3',
      ],
  python_requires='>=3',
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Science/Research',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT',   # Again, pick a license
  ],
)