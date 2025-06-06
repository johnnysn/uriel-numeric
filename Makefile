# Root Makefile

# Default target
all: numeric tools cellsolver monodomainsolver

numeric:
	@echo "Building uriel-numeric..."
	@cmake -S uriel-numeric -B uriel-numeric
	@cmake --build uriel-numeric

tools:
	@echo "Building uriel-tools..."
	@cmake -S uriel-tools -B uriel-tools
	@cmake --build uriel-tools

cellsolver:
	@echo "Building cardiac-cell-solver..."
	@cmake -S cardiac-cell-solver -B cardiac-cell-solver
	@cmake --build cardiac-cell-solver

monodomainsolver:
	@echo "Building cardiac-monodomain-fd-solver..."
	@cmake -S cardiac-monodomain-fd-solver -B cardiac-monodomain-fd-solver
	@cmake --build cardiac-monodomain-fd-solver

clean:
	@echo "Cleaning all ..."
	@rm -rf uriel-numeric
	@rm -rf uriel-tools
	@rm -rf cardiac-cell-solver
	@rm -rf cardiac-monodomain-fd-solver

.PHONY: all numeric tools cellsolver monodomainsolver clean