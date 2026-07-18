# generate_index.py needs the pixi environment's Python (3.10+, for the union
# type hints); the pixi.toml manifest lives in the project root, one folder up.
# --manifest-path selects that environment while keeping the working directory
# here in view/, so --output and the */*.html scan stay relative to this folder.
PIXI := pixi run --manifest-path ../pixi.toml

REPORTS := $(wildcard */*.html)

.PHONY: all force clean

all: index.html

index.html: generate_index.py $(REPORTS)
	$(PIXI) python generate_index.py --output $@

# Regenerate every per-folder report.html from the report payloads (overwriting
# existing ones), then rebuild the top-level index.
force: generate_index.py $(REPORTS)
	$(PIXI) python generate_index.py --force --output index.html

clean:
	rm -f index.html
