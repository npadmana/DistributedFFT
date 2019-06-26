Mason.toml: Mason.toml.template
	envsubst < $< > $@

target/example/plane_v2: example/Comm/plane_v2.chpl
	chpl -o $@ $< --fast -schpl_serializeSlices -suseBulkTransfer
