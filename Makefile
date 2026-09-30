JULIA ?= julia

.PHONY: deploy

deploy:
	$(JULIA) --project=. -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=. app/GadgetEditor.jl
