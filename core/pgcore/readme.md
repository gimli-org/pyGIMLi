Core extension for [pyGIMLi](www.pygimli.org)

This package is supposed to be the platform dependent part of the 'any-platform' `pygimli` package.
If you install `pygimli` via pip or conda, a suitable `pgcore` package will be installed as dependency. 

If we post small bugfixes or repair some issues with fast git commits, you can use the `pgcore` package as fallback core until we publish a new release.
`pygimli` will try to import manual builded cores or search for global installed cores from `PYTHONPATH`.
