EXAMPLE_SRC=$(shell find example -name *.chpl)

EXAMPLES=target/example/fftw-mpi-benchmark \
	target/example/plane_v2 \
	target/example/FFTW/fftw-plan-timings \
	target/example/Dist/time_warmup



Mason.toml: Mason.toml.template
	envsubst < $< > $@

target/example/plane_v2: example/Comm/plane_v2.chpl
	chpl -o $@ $< --fast -schpl_serializeSlices -suseBulkTransfer

target/example/fftw-mpi-benchmark: example/MPI/fftw-mpi-benchmark.chpl src/DistributedFFT.chpl
	chpl -o $@ $< --fast ${MPI_CHPL_FLAGS} ${CHPL_WARN_FLAGS} -lfftw3_mpi -lfftw3_threads -lfftw3 -Msrc

target/example/FFTW/%: example/FFTW/%.chpl
	mkdir -p target/example/FFTW
	chpl -o $@ $< --fast ${MPI_CHPL_FLAGS} ${CHPL_WARN_FLAGS} -lfftw3_threads -lfftw3 --local

target/example/Dist/%: example/Dist/%.chpl src/DistributedFFT.chpl
	mkdir -p target/example/Dist
	chpl -o $@ $< --fast ${CHPL_WARN_FLAGS} -lfftw3_threads -lfftw3 -Msrc


.PHONY: examples target
examples: target ${EXAMPLES}

target:
	mkdir -p target/example

