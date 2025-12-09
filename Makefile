MODULE_TOPDIR = ../..

PGM = i.hyper.continuum

include $(MODULE_TOPDIR)/include/Make/Script.make
include $(MODULE_TOPDIR)/include/Make/Html.make

default: script html
