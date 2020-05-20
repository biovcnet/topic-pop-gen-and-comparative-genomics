import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
import pyfaidx as px
def degenerancy(data,codonDict):
	
	#DEGENERANCY DICTIONARIES
	standardDict = {
	'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '022', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	VertebrateMitochondrial = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '004', 'CGG': '004', 'ATT': '002', 'ATC': '002', 'ATA': '002', 'ATG': '002', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '002', 'AGG': '002', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	YeastMitochondrial = {

		'TTT': '002', 'TTC': '002', 'TTA': '002', 'TTG': '002', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '004', 'CTG': '004', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '002', 'ATC': '002', 'ATA': '002', 'ATG': '002', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	MoldMitochondrialProtozoanMitochondrialCoelenterateMitochondrialMycoplasmaSpiroplasma = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	InvertebrateMitochondrial = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '004', 'CGG': '004', 'ATT': '002', 'ATC': '002', 'ATA': '002', 'ATG': '002', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '004', 'AGC': '004', 'AGA': '004', 'AGG': '004', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	CiliateNuclearDasycladaceanNuclearHexamitaNuclear = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '202', 'TAG': '202', 'TGT': '002', 'TGC': '002', 'TGA': '000', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '202', 'CAG': '202', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	EchinodermMitochondrialFlatwormMitochondrial = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '004', 'CGG': '004', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '003', 'AAC': '003', 'AAA': '003', 'AAG': '000', 'AGT': '004', 'AGC': '004', 'AGA': '004', 'AGG': '004', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	EuplotidNuclear = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '003', 'TGC': '003', 'TGA': '003', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	BacterialArchaealandPlantPlastid = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '022', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	AlternativeYeastNuclear = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '002', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '022', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000', 'CTT': '003', 'CTC': '003', 'CTA': '203', 'CTG': '000', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	AlternativeFlatwormMitochondrial = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '003', 'TAC': '003', 'TAA': '003', 'TAG': '000', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '004', 'CGG': '004', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '003', 'AAC': '003', 'AAA': '003', 'AAG': '000', 'AGT': '004', 'AGC': '004', 'AGA': '004', 'AGG': '004', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	BlepharismaMacronuclear = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '020', 'TAG': '200', 'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '202', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	ChlorophyceanMitochondrial = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '222', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '020', 'TAG': '020', 'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	TrematodeMitochondrial = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '004', 'CGG': '004', 'ATT': '002', 'ATC': '002', 'ATA': '002', 'ATG': '002', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '003', 'AAC': '003', 'AAA': '003', 'AAG': '000', 'AGT': '004', 'AGC': '004', 'AGA': '004', 'AGG': '004', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	ScenedesmusobliquusMitochondrial = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '222', 'TCT': '003', 'TCC': '003', 'TCA': '030', 'TCG': '003', 'TAT': '002', 'TAC': '002', 'TAA': '030', 'TAG': '020', 'TGT': '002', 'TGC': '002', 'TGA': '030', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	ThraustochytriumMitochondrial = {

		'TTT': '002', 'TTC': '002', 'TTA': '030', 'TTG': '200', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '032', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '030', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '004', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	PterobranchiaMitochondrial = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '004', 'CGG': '004', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '022', 'AGT': '003', 'AGC': '003', 'AGA': '003', 'AGG': '020', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	CandidateDivisionSR1andGracilibacteria = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '200', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '204', 'GGG': '004'}
	PachysolentannophilusNuclear = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '002', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '022', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000', 'CTT': '003', 'CTC': '003', 'CTA': '203', 'CTG': '000', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	KaryorelictNuclear = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '202', 'TAG': '202', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '202', 'CAG': '202', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	CondylostomaNuclear = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '202', 'TAG': '202', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '202', 'CAG': '202', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	MesodiniumNuclear = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '004', 'TAC': '004', 'TAA': '004', 'TAG': '004', 'TGT': '002', 'TGC': '002', 'TGA': '000', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	PeritrichNuclear = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '202', 'TAG': '202', 'TGT': '002', 'TGC': '002', 'TGA': '000', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '202', 'GAG': '202', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	BlastocrithidiaNuclear = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '202', 'TAG': '202', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '202', 'GAG': '202', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
	BalanophoraceaePlastid = {

		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '000', 'TAG': '000', 'TGT': '002', 'TGC': '002', 'TGA': '000', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}

	if(codonDict == 'standard'):
		degenerateCodonTable = standardDict
	elif(codonDict == 'VertebrateMitochondrial'):
		degenerateCodonTable = VertebrateMitochondrial
	elif(codonDict == 'YeastMitochondrial'):
		degenerateCodonTable = YeastMitochondrial
	elif(codonDict == 'MoldMitochondrialProtozoanMitochondrialCoelenterateMitochondrialMycoplasmaSpiroplasma'):
		degenerateCodonTable = MoldMitochondrialProtozoanMitochondrialCoelenterateMitochondrialMycoplasmaSpiroplasma
	elif(codonDict == 'InvertebrateMitochondrial'):
		degenerateCodonTable = InvertebrateMitochondrial
	elif(codonDict == 'CiliateNuclearDasycladaceanNuclearHexamitaNuclear'):
		degenerateCodonTable = CiliateNuclearDasycladaceanNuclearHexamitaNuclear
	elif(codonDict == 'EchinodermMitochondrialFlatwormMitochondrial'):
		degenerateCodonTable = EchinodermMitochondrialFlatwormMitochondrial
	elif(codonDict == 'EuplotidNuclear'):
		degenerateCodonTable = EuplotidNuclear
	elif(codonDict == 'BacterialArchaealandPlantPlastid'):
		degenerateCodonTable = BacterialArchaealandPlantPlastid
	elif(codonDict == 'AlternativeYeastNuclear'):
		degenerateCodonTable = AlternativeYeastNuclear
	elif(codonDict == 'AlternativeFlatwormMitochondrial'):
		degenerateCodonTable = AlternativeFlatwormMitochondrial
	elif(codonDict == 'BlepharismaMacronuclear'):
		degenerateCodonTable = BlepharismaMacronuclear
	elif(codonDict == 'ChlorophyceanMitochondrial'):
		degenerateCodonTable = ChlorophyceanMitochondrial
	elif(codonDict == 'TrematodeMitochondrial'):
		degenerateCodonTable = TrematodeMitochondrial
	elif(codonDict == 'ScenedesmusobliquusMitochondrial'):
		degenerateCodonTable = ScenedesmusobliquusMitochondrial	
	elif(codonDict == 'ThraustochytriumMitochondrial'):
		degenerateCodonTable = ThraustochytriumMitochondrial
	elif(codonDict == 'PterobranchiaMitochondrial'):
		degenerateCodonTable = PterobranchiaMitochondrial	
	elif(codonDict == 'CandidateDivisionSR1andGracilibacteria'):
		degenerateCodonTable = CandidateDivisionSR1andGracilibacteria
	elif(codonDict == 'PachysolentannophilusNuclear'):
		degenerateCodonTable = PachysolentannophilusNuclear
	elif(codonDict == 'KaryorelictNuclear'):
		degenerateCodonTable = KaryorelictNuclear
	elif(codonDict == 'CondylostomaNuclear'):
		degenerateCodonTable = CondylostomaNuclear
	elif(codonDict == 'MesodiniumNuclear'):
		degenerateCodonTable = MesodiniumNuclear
	elif(codonDict == 'PeritrichNuclear'):
		degenerateCodonTable = PeritrichNuclear
	elif(codonDict == 'BlastocrithidiaNuclear'):
		degenerateCodonTable = BlastocrithidiaNuclear
	elif(codonDict == 'BalanophoraceaePlastid'):
		degenerateCodonTable = BalanophoraceaePlastid



	degenerancy = ''
	for i in range(0, len(data),3):
		codon = data[i:i+3]
		if('N' in codon or '-' in codon):
			degenerancy += codon
		else:
			degenerancy += degenerateCodonTable[codon]

	return(degenerancy)
def sequencesToMatrix(multiFasta,split=None):

	# Extract samples from fastas
	samples = list(multiFasta.keys())
	
	if(split is None):
		seqLen = len(multiFasta[samples[0]][:].seq)
		
		if((seqLen % 3) != 0):
			print('cdsLength')
			sys.exit('cdsLength')

		# Create empty array with ndimesions equal to multi-Fasta lines and length
		matrix = np.empty([len(samples),len(multiFasta[samples[0]][:].seq)],dtype='str')
		
		# List to append indexes if whole sequence at any population is len(seq) * 'N'
		deleteIndex = list()
		
		# Iter fasta to add sequence to matrix
		for i in range(1,len(samples),1):
			# Extract each sample sequence
			tmp = multiFasta[samples[i]][:].seq
			if(len(tmp) != seqLen):
				print('errorAlign')
				sys.exit('errorAlign')

			# if('N' in tmp):
			if('N'*len(tmp) == tmp):
				deleteIndex.append(i)
			else:
				matrix[i] = list(tmp)
		
		# Delete lines
		matrix = np.delete(matrix,deleteIndex,0)
		
		degenCode = degenerancy(multiFasta[samples[0]][:].seq,args.codonTable)
		# Put degenerancy in first ndarray element
		matrix[0] = list(degenCode)
		# NEED TO SOLVE THIS. C-contigous change web subset, need true in order to inter properly witih nditer

		matrix = np.asarray(matrix[:,(matrix[0]=='0') | (matrix[0]=='4')],order='C')

		return(matrix)

	else:
		# Create empty array with ndimesions equal to multi-Fasta lines and length
		splitC = [split[0],split[1]]
		matrix = np.empty([len(samples),(splitC[1]-splitC[0]+1)],dtype='str')
		
		# List to append indexes if whole sequence at any population is len(seq) * 'N'
		deleteIndex = list()
		
		# Iter fasta to add sequence to matrix
		for i in range(1,len(samples),1):
			# Extract each sample sequence
			tmp = multiFasta.get_spliced_seq(samples[i], [splitC]).seq
			if(tmp == ('N' * len(tmp))):
				deleteIndex.append(i)
			else:
				matrix[i] = list(tmp)
		
		# Delete lines
		matrix = np.delete(matrix,deleteIndex,0)
		
		degenCode = degenerancy(multiFasta.get_spliced_seq(samples[0], [splitC]).seq.upper(),args.codonTable)
		# Put degenerancy in first ndarray element
		matrix[0] = list(degenCode)
		# NEED TO SOLVE THIS. C-contigous change web subset, need true in order to inter properly witih nditer

		matrix = np.asarray(matrix[:,(matrix[0]=='0') | (matrix[0]=='4')],order='C')  

		return(matrix)
def uSfsFromFasta(sequenceMatrix):
	output = list()
	for x in np.nditer(sequenceMatrix, order='F',flags=['external_loop']): 
		degen = x[0]
		AA = x[-1]

		# Undefined Ancestra Allele. Try to clean out of the loop
		if(AA == 'N' or AA == '-'):
			next
		elif('N' in x[1:-1] or '-' in x[1:-1]):
			next
		# Monomorphic sites. Try to clean out of the loop
		elif(np.unique(x[1:][np.where(x[1:]!='N')]).shape[0] == 1):
			next
		else:
			pol = x[1:-1]
			if(degen == '4'):
				functionalClass = '4fold'
			else:
				functionalClass = '0fold'

			# Check if pol != AA and monomorphic
			if((np.unique(pol).shape[0] == 1) and (np.unique(pol)[0] != AA)):
				div = 1; AF = 0
				tmp = [AF,div,functionalClass]
				output.append(tmp)
			else:
				AN = x[1:-1].shape[0]
				AC = pd.DataFrame(data=np.unique(x[1:-1], return_counts=True)[1],index=np.unique(x[1:-1], return_counts=True)[0])
				div = 0
				if(AA not in AC.index):
					next
				else:
					AC = AC[AC.index!=AA]
					if(len(AC) == 0):
						next
					else:
						AF = AC.iloc[0]/AN
						AF = AF.iloc[0]
				tmp = [AF,div,functionalClass]
				output.append(tmp)
	return(output)
def formatSfs(sequenceMatrix,rawSfsOutput,dafFile,divFile,path,append=True):

	df = pd.DataFrame(rawSfsOutput)
	df['id'] = 'uploaded'
	df.columns = ['derivedAlleleFrequency','d','functionalClass','id']

	# Extract divergence data
	div = df[['id','functionalClass','d']]
	div = div[div['d']!=0]
	div = div.groupby(['id','functionalClass'])['d'].count().reset_index()
	div = div.pivot_table(index=['id'],columns=['functionalClass'],values='d').reset_index()
	try:
		div = div[['0fold','4fold']]
	except:
		if('4fold' in div.columns):
			div = div[['4fold']]
			div['0fold'] = 0
		elif('0fold' in div.columns):
			div = div[['0fold']]
			div['4fold'] = 0

	div['mi'] =  sequenceMatrix[0][np.where(sequenceMatrix[0]=='0')].shape[0]
	div['m0'] =  sequenceMatrix[0][np.where(sequenceMatrix[0]=='4')].shape[0]
	div = div.rename(columns={'0fold':'Di','4fold':'D0','mi':'mi','m0':'m0'})
	# div = div.pivot_table(index=['functionalClass'],columns=['functionalClass'],values='div').reset_index()

	# Create SFS pd.DataFrame by functionClass and 20 frequency bin
	daf = df[df['d']!=1][['derivedAlleleFrequency','functionalClass','id']]
	bins = np.arange(0.025,1.05,0.05)
	labels = np.arange(0.025,1.0,0.05).tolist()
	daf['categories'] = pd.cut(daf['derivedAlleleFrequency'],bins=bins,labels=labels)
	daf = daf.groupby(['functionalClass','id','categories']).count().reset_index()
	sfs = pd.DataFrame({'daf':daf['categories'].unique(),'P0':daf[daf['functionalClass']=='4fold']['derivedAlleleFrequency'].reset_index(drop=True),'Pi':daf[daf['functionalClass']=='0fold']['derivedAlleleFrequency'].reset_index(drop=True)})

	sfs = sfs[['daf','P0','Pi']]
	sfs['P0'] = sfs['P0'].fillna(0)
	sfs['Pi'] = sfs['Pi'].fillna(0)
	sfs['daf'] = sfs['daf'].apply(lambda x: round(x,3))

	if(append is True):
		sfs.to_csv(path + dafFile,sep='\t',header=True,index=False,mode='a')
		div.to_csv(path + divFile,sep='\t',header=True,index=False,mode='a')
	else:
		sfs.to_csv(path + dafFile,sep='\t',header=True,index=False)
		div.to_csv(path + divFile,sep='\t',header=True,index=False)

#################BGD GROUP######################
########UNIVERSITAT AUTONOMA DE BARCELONA#######
# We highly recommend to read the code and test manual files before run large analysis
if __name__ == '__main__':
	'''Parse arguments and show the required inputs if only name is given to command line'''
	parser = argparse.ArgumentParser(description='Estimate binned DAF and diverngece from multiFASTA alignments. The expected input is the same as described in https://doi.org/10.1093/nar/gkz372')
	# Required arguments
	parser.add_argument('--multiFasta', type = str, required = True, help = 'Raw data to estimate SFS and divergence, including reference and outgroup.')
	parser.add_argument('--daf', type = str, required = True, help = 'Name to DAF file')
	parser.add_argument('--div', type = str, required = True, help = 'Name to divergence file')
	parser.add_argument('--codonTable', type = str, required = True, choices=['standard','VertebrateMitochondrial','YeastMitochondrial','MoldMitochondrialProtozoanMitochondrialCoelenterateMitochondrialMycoplasmaSpiroplasma','InvertebrateMitochondrial','CiliateNuclearDasycladaceanNuclearHexamitaNuclear','EchinodermMitochondrialFlatwormMitochondrial','EuplotidNuclear','BacterialArchaealandPlantPlastid','AlternativeYeastNuclear','AlternativeFlatwormMitochondrial','BlepharismaMacronuclear','ChlorophyceanMitochondrial','TrematodeMitochondrial','ScenedesmusobliquusMitochondrial','ThraustochytriumMitochondrial','PterobranchiaMitochondrial','CandidateDivisionSR1andGracilibacteria','PachysolentannophilusNuclear','KaryorelictNuclear','CondylostomaNuclear','MesodiniumNuclear','PeritrichNuclear','BlastocrithidiaNuclear','BalanophoraceaePlastid'],help = 'Codon degeneracy to use')
	parser.add_argument('--startCoordinates', type = str, required = False, help = 'List of coordinates if CDS are merge in one large CDS')

	args = parser.parse_args()
	start = time.time()
	pwd = os.getcwd() + '/'

	# Open multi-Fasta
	file = px.Fasta(args.multiFasta,duplicate_action='first',sequence_always_upper=True,read_long_names=True)

	# Create ndarray with sequences
	multiFastaMatrix = sequencesToMatrix(file)

	# Check if there is more than 2 individuals to extract polymorphism
	if(multiFastaMatrix.shape[0] < 4):
		print('numberOfLines')
		sys.exit('numberOfLines')

	# Estimating SFS
	rawSfs = uSfsFromFasta(multiFastaMatrix)

	# Formating SFS
	formatSfs(multiFastaMatrix,rawSfs,args.daf,args.div,pwd)





