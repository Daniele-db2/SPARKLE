#!/usr/bin/python env

#python 3 standard library

import argparse
from argparse import HelpFormatter
import sys

def main():

	parser = argparse.ArgumentParser(prog='SPARKLE', description='''Proof of concept: alignment of FASTQ sequences to a refernece genome using clustered SPARK''', epilog='''This program was developed by Daniele Bellini (https://github.com/Daniele-db2)''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='index,align')

	## index ##

	parser_index = subparsers.add_parser('index', help='Create hash table with user-defined k-mer size for a reference genome')

	required = parser_index.add_argument_group('Required I/O arguments')

	required.add_argument('-g', '--genome', help='reference genome', metavar='FASTA', required=True)

	additional = parser_index.add_argument_group('Additional parameters')

	additional.add_argument('-k','--kmer', help='k-mer size for indexing [10]', default=10, type=int)

	parser_index.set_defaults(func=run_subtool)
	
	## align ##

	parser_align = subparsers.add_parser('align', help='Align sequences in FASTQ format to an indexed reference genome')

	required = parser_align.add_argument_group('Required I/O arguments')

	required.add_argument('-g', '--genome', help='reference genome', metavar='FASTA', required=True)
	required.add_argument('-r', '--reads', help='sequencing reads', metavar='FASTQ', required=True)
	required.add_argument('-o', '--output', help='output compressed alignment in text format', metavar='TXT.GZ', required=True)

	additional = parser_align.add_argument_group('Additional parameters')

	additional.add_argument('-k','--kmer', help='k-mer size for indexing [10]', default=10, type=int)
	additional.add_argument('--match', help='match reward (global/local alignment) [5]', default=5, type=int)
	additional.add_argument('--mismatch', help='mismatch penalty (global/local alignment) [-4]', default=-4, type=int)
	additional.add_argument('--gapopen', help='gap opening penalty (global/local alignment) [-8]', default=-8, type=int)
	additional.add_argument('--gapextend', help='gap extending penalty (global alignment) [-6]', default=-6, type=int)
	additional.add_argument('--alignment', help='which alignment to use. Choose between global, local and guess [guess]', choices=['global', 'local', 'guess'], default='guess', type=str)
	additional.add_argument('--distance', help='minium distance for chaining [50]', type=int, default=50)

	parser_align.set_defaults(func=run_subtool)


	#print help if no subcommand nor --help provided
	
	if len(sys.argv)==1:
		
		parser.print_help(sys.stderr)
		sys.exit(1)

	#case-insensitive submodules
	
	if sys.argv[1].lower() == 'index':

		sys.argv[1] = 'index'

	elif sys.argv[1].lower() == 'align':

		sys.argv[1] = 'align'

	args = parser.parse_args()
	args.func(parser, args)



class CustomFormat(HelpFormatter):


	def _format_action_invocation(self, action):

		if not action.option_strings:

			default = self._get_default_metavar_for_positional(action)
			metavar, = self._metavar_formatter(action, default)(1)
			
			return metavar

		else:

			parts = []

			if action.nargs == 0:

				parts.extend(action.option_strings)

			else:

				default = self._get_default_metavar_for_optional(action)
				args_string = self._format_args(action, default)
				
				for option_string in action.option_strings:

					parts.append(option_string)

				return '%s %s' % (', '.join(parts), args_string)

			return ', '.join(parts)

	def _get_default_metavar_for_optional(self, action):

		return action.dest.upper()



def run_subtool(parser, args):


	if args.command == 'index': #index

		from .index import index as submodule

	elif args.command == 'align': #align

		from .align import align as submodule

	else:

		parser.print_help()

	submodule.run(parser,args)


if __name__ =='__main__':


	main()