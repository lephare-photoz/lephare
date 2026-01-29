SPHINXOPTS    ?= -T -E
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = ../_readthedocs/html
NOTEBOOK_DIR  = notebooks
NOTEBOOKS     := $(wildcard $(NOTEBOOK_DIR)/*.ipynb)

.PHONY: help clean html no-nb notebooks force-notebooks

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Reconstruit toute la documentation (y compris les notebooks modifiés)
html: $(NOTEBOOKS)
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Reconstruit sans les notebooks
no-nb:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) -D exclude_patterns="$(NOTEBOOK_DIR)/*" $(O)

# Reconstruit uniquement les notebooks modifiés
notebooks: $(NOTEBOOKS)

# Règle pour chaque notebook : reconstruction si modifié
$(NOTEBOOK_DIR)/%.ipynb:
	@echo "Reconstruction du notebook : $<"
	@jupyter nbconvert --to notebook --inplace --execute $<

# Force la reconstruction de tous les notebooks
force-notebooks:
	@for notebook in $(NOTEBOOKS); do \
		jupyter nbconvert --to notebook --inplace --execute "$$notebook"; \
	done

clean:
	rm -rf "$(BUILDDIR)" _build/doctrees
