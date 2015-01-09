SUBDIRS = pipeline/qpid pipeline/brahle_assembly pipeline/msa

.PHONY: components $(SUBDIRS)

components: clean prepare_bin $(SUBDIRS) link_bin

prepare_bin:
	@test -d bin || mkdir bin

link_bin:
	ln -s ../pipeline/qpid/bin/overlap bin/croler_overlap
	ln -s ../pipeline/brahle_assembly/bin/main_layout bin/croler_layout
	ln -s ../pipeline/msa/bin/msa bin/croler_consensus

clean:
	@if [ -d bin ]; then rm -r bin; fi;
	@make -C pipeline/qpid clean
	@make -C pipeline/brahle_assembly clean
	@make -C pipeline/msa clean

$(SUBDIRS):
	$(MAKE) -C $@
