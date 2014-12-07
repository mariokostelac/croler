SUBDIRS = pipeline/qpid pipeline/brahle_assembly pipeline/msa

.PHONY: components $(SUBDIRS)

components: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@
