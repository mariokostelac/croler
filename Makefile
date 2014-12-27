SUBDIRS = pipeline/qpid pipeline/brahle_assembly pipeline/msa

.PHONY: components $(SUBDIRS)

components: prepare_bin $(SUBDIRS) copy_bin

prepare_bin:
	@test -d bin || mkdir bin

copy_bin:
	cp pipeline/qpid/bin/overlap bin/overlap
	cp pipeline/brahle_assembly/bin/main_layout bin/layout
	cp pipeline/msa/bin/msa bin/consensus

$(SUBDIRS):
	$(MAKE) -C $@
