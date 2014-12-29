SUBDIRS = pipeline/qpid pipeline/brahle_assembly pipeline/msa

.PHONY: components $(SUBDIRS)

components: clean prepare_bin $(SUBDIRS) link_bin

prepare_bin:
	@test -d bin || mkdir bin

link_bin:
	ln -s ../pipeline/qpid/bin/overlap bin/overlap
	ln -s ../pipeline/brahle_assembly/bin/main_layout bin/layout
	ln -s ../pipeline/msa/bin/msa bin/consensus

clean:
	@if [ -d bin ]; then rm -r bin; fi;

$(SUBDIRS):
	$(MAKE) -C $@
