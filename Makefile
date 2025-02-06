# Root Makefile

# Default target (if no specific target is given, this will run)
all: numeric tools cellsolver

numeric:
	@echo "Building uriel-numeric..."
	$(MAKE) -C uriel-numeric

tools:
	@echo "Building uriel-tools..."
	$(MAKE) -C uriel-tools

cellsolver:
	@echo "Building cardiac-cell-solver..."
	$(MAKE) -C cardiac-cell-solver

monodomainsolver:
	@echo "Building cardiac-monodomain-fd-solver..."
	$(MAKE) -C cardiac-monodomain-fd-solver

clean:
	@echo "Cleaning up uriel-numeric..."
	$(MAKE) clean -C uriel-numeric
	@echo "Cleaning up uriel-tools..."
	$(MAKE) clean -C uriel-tools
	@echo "Cleaning up cardiac-cell-solver..."
	$(MAKE) clean -C cardiac-cell-solver
	@echo "Cleaning up cardiac-monodomain-fd-solver..."
	$(MAKE) clean -C cardiac-monodomain-fd-solver

.PHONY: all numeric tools clean