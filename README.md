# cornish
A pure Python interface to the excellent [Starlink AST library](http://starlink.eao.hawaii.edu/starlink/AST). The C library handles WCS systems for astronomical data, including conversion between systems and plotting (and a lot more!).

There is an existing [Python interface](https://pypi.python.org/pypi/starlink-pyast/3.9.0) actively supported by the authors (David Berry and Tim Jenness), but it effectively only exposes C functions in Python. You must be very familiar with the C version of the library to make use of this (I recommend keeping the [documentation](http://starlink.eao.hawaii.edu/devdocs/sun211.htx/sun211.html) handy!).

Cornish is a new Python package written to be a fully Pythonic interface to the library to make it more accessible, self-documenting, and easy to use.

This is a very early work in progress; there is no working code here yet. Contributions welcome.

