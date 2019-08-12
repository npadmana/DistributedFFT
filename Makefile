EXAMPLE_SRC=$(shell find example -name *.chpl)

EXAMPLES=target/example/fftw-mpi-benchmark \
	target/example/plane_v2 \
	target/example/FFTW/fftw-plan-timings \
	target/example/Dist/time_warmup \
	target/example/Dist/test_simple \
	target/example/Dist/time_simple \


# The following are experimental flags for better performance
PERF_FLAGS=-schpl_serializeSlices -suseBulkTransfer


Mason.toml: Mason.toml.template
	envsubst < $< > $@

target/example/plane_v2: example/Comm/plane_v2.chpl
	chpl -o $@ $< --fast ${PERF_FLAGS}

target/example/fftw-mpi-benchmark: example/MPI/fftw-mpi-benchmark.chpl src/DistributedFFT.chpl
	chpl -o $@ $< --fast ${MPI_CHPL_FLAGS} ${CHPL_WARN_FLAGS} -lfftw3_mpi -lfftw3_threads -lfftw3 -Msrc

target/example/FFTW/%: example/FFTW/%.chpl
	mkdir -p target/example/FFTW
	chpl -o $@ $< --fast ${MPI_CHPL_FLAGS} ${CHPL_WARN_FLAGS} -lfftw3_threads -lfftw3 --local

target/example/Dist/%: example/Dist/%.chpl src/DistributedFFT.chpl
	mkdir -p target/example/Dist
	chpl -o $@ $< --fast ${CHPL_WARN_FLAGS} -lfftw3_threads -lfftw3 -Msrc ${PERF_FLAGS}

target/example/NPB-FT/%: example/NPB-FT/%.chpl src/DistributedFFT.chpl
	mkdir -p target/example/NPB-FT
	chpl -o $@ $< --fast ${CHPL_WARN_FLAGS} -lfftw3 -Msrc ${PERF_FLAGS}

ft: target/example/NPB-FT/ft target/example/NPB-FT/ft_transposed

target/example/R2R/%: example/R2R/%.chpl src/DistributedFFT.chpl
	mkdir -p target/example/R2R
	chpl -o $@ $< ${CHPL_WARN_FLAGS} -lfftw3 -Msrc ${PERF_FLAGS}

.PHONY: examples target
examples: target ${EXAMPLES}

target:
	mkdir -p target/example

