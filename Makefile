Mason.toml: Mason.toml.template
	envsubst < $< > $@

target/example/plane_v2: example/Comm/plane_v2.chpl
	chpl -o $@ $< --fast -schpl_serializeSlices -suseBulkTransfer

target/example/fftw-mpi-benchmark: example/MPI/fftw-mpi-benchmark.chpl src/DistributedFFT.chpl
	chpl -o $@ $< --fast ${MPI_CHPL_FLAGS} ${CHPL_WARN_FLAGS} -lfftw3_mpi -lfftw3_threads -lfftw3 -Msrc/


EXAMPLES=target/example/fftw-mpi-benchmark \
	target/example/plane_v2


.PHONY: examples
examples: ${EXAMPLES}
