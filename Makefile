Mason.toml: Mason.toml.template
	envsubst < $< > $@
