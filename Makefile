NAUTY_VERSION := 2.9.3
NAUTY_DIR     := nauty-$(NAUTY_VERSION)
NAUTY_TAR     := nauty$(subst .,,$(NAUTY_VERSION)).tar.gz
NAUTY_URL     := https://pallini.di.uniroma1.it/$(NAUTY_TAR)
PREFIX        ?= $(CURDIR)/deps

.PHONY: setup clean check

setup: $(PREFIX)/bin/shortg
	@echo "shortg installed at $(PREFIX)/bin/shortg"
	@echo "Add to PATH:  export PATH=$(PREFIX)/bin:\$$PATH"

check:
	@if command -v shortg >/dev/null 2>&1; then \
		echo "shortg found: $$(command -v shortg)"; \
		shortg --version 2>&1 | head -1; \
	else \
		echo "shortg not found in PATH. Run 'make setup' to build from source."; \
	fi

$(PREFIX)/bin/shortg: $(NAUTY_DIR)/shortg
	mkdir -p $(PREFIX)/bin
	cp $(NAUTY_DIR)/shortg $(PREFIX)/bin/shortg

$(NAUTY_DIR)/shortg: $(NAUTY_DIR)/configure
	cd $(NAUTY_DIR) && ./configure --quiet && $(MAKE) shortg

$(NAUTY_DIR)/configure: $(NAUTY_TAR)
	tar xzf $(NAUTY_TAR)

$(NAUTY_TAR):
	curl -fLO $(NAUTY_URL)

clean:
	rm -rf $(NAUTY_DIR) $(NAUTY_TAR)
