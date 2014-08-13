#!/bin/bash

runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_500_97 -D MIN_LEN=500 -D ALN_WIGGLE=15 -D MIN_IDENT=97 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_500_95 -D MIN_LEN=500 -D ALN_WIGGLE=15 -D MIN_IDENT=95 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_500_93 -D MIN_LEN=500 -D ALN_WIGGLE=15 -D MIN_IDENT=93 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_500_91 -D MIN_LEN=500 -D ALN_WIGGLE=15 -D MIN_IDENT=91 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_1000_97 -D MIN_LEN=1000 -D ALN_WIGGLE=15 -D MIN_IDENT=97 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_1000_95 -D MIN_LEN=1000 -D ALN_WIGGLE=15 -D MIN_IDENT=95 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_1000_93 -D MIN_LEN=1000 -D ALN_WIGGLE=15 -D MIN_IDENT=93 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_1000_91 -D MIN_LEN=1000 -D ALN_WIGGLE=15 -D MIN_IDENT=91 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_2000_97 -D MIN_LEN=2000 -D ALN_WIGGLE=15 -D MIN_IDENT=97 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_2000_95 -D MIN_LEN=2000 -D ALN_WIGGLE=15 -D MIN_IDENT=95 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_2000_93 -D MIN_LEN=2000 -D ALN_WIGGLE=15 -D MIN_IDENT=93 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_2000_91 -D MIN_LEN=2000 -D ALN_WIGGLE=15 -D MIN_IDENT=91 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_3000_97 -D MIN_LEN=3000 -D ALN_WIGGLE=15 -D MIN_IDENT=97 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_3000_95 -D MIN_LEN=3000 -D ALN_WIGGLE=15 -D MIN_IDENT=95 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_3000_93 -D MIN_LEN=3000 -D ALN_WIGGLE=15 -D MIN_IDENT=93 -D FASTA_EXP=1 preassembled_reads.fasta
runAmos -C ${PBX_PROTOCOL_TEMPLATES_PATH}/MinimoNucmer.acf -D PBX_MUMMER_PATH=${PBX_MUMMER_PATH} -D QUAL_IN=preassembled_reads.qual -D OUT_PREFIX=tale_assembly_3000_91 -D MIN_LEN=3000 -D ALN_WIGGLE=15 -D MIN_IDENT=91 -D FASTA_EXP=1 preassembled_reads.fasta