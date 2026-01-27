# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = ./
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Used to auto-update website in local build. Requires sphinx-autobuild.
auto:
	sphinx-autobuild --delay 2 --ignore "pygimliapi/*" --ignore "pygimliapi/_generated/*" --ignore "*gallery.html" --ignore "*.examples" --ignore "*.ipynb" --ignore "*.zip" --ignore "*.pickle" . _build/html -D plot_gallery=0

clean:
	rm -rf $(BUILDDIR)
	rm -rf _examples_auto
	rm -rf _tutorials_auto
	rm -rf pygimliapi/_generated
	rm -rf examples/modelling/*.geo*
	rm -rf examples/modelling/*.msh*
	rm -f _templates/gallery.html
	rm -rf jupyter_execute
	rm -rf .jupyter_cache

doc: html

api:
	mkdir -p pygimliapi/_generated
	cp -f pygimliapi/viewer_overwrite pygimliapi/_generated/pygimli.viewer.rst
	cp -f pygimliapi/physics_overwrite pygimliapi/_generated/pygimli.physics.rst

html-nogallery: api
	$(SPHINXBUILD) -D plot_gallery=0 -b dirhtml "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished without gallery. The HTML pages are in $(BUILDDIR)."

html: api
	$(SPHINXBUILD) -v -D abort_on_example_error=True -b html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)."

serve:
	python -m http.server 11000 --directory _build
