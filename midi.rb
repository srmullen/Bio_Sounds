#!usr/local/bin/ruby
#
# Programmer: Sean Mullen
# 
# == Description
#
# MIDI Format class
#
#=begin
require 'bio/db'
require 'bio/sequence'
require 'bio/sequence/dblink'
require 'bio/db/fasta/defline'
require 'bio/db/fasta'
#=end
module Bio
	class MidiFormat < FastaFormat

		attr_reader :aa_first_pos, :aa_second_pos, :aa_third_pos, :midi_table, :midi_data, :tempo, :track_data
		attr_writer :tempo

		def initialize(str)
			super(str)

			@tempo = 77 # will be up to the user.

			@midi_table = Hash.new
			@na_hash = Hash.new
			@aa_hash = Hash.new

			# Create an array of arrays to store the result of
			# compute_midi_data, and preallocate them to the length
			# of the sequence.
			@midi_notes = Array.new
			@midi_notes.push(Array.new)
			@midi_notes.push(Array.new)
			@midi_notes.push(Array.new)

			# Create an array to hold the completed track data
			@track_data = Array.new

			create_aminos
			create_midi_table
			compute_midi_notes
			create_track_data(0)
			create_track_data(1)
			create_track_data(2)
		end
		
		# Method: create_aminos
		# => Translates the NA sequence into an AA sequence for
		# use in determining the midi number.
		def create_aminos
			nas = naseq
			@aa_first_pos = nas.translate.to_s.strip
			@aa_second_pos = nas.translate(2).to_s.strip
			@aa_third_pos = nas.translate(3).to_s.strip
		end

		# Method: create_midi_hash
		# => Creates a hash of hashes so a midi number can be
		# indexed by its associated nucleic acid and amino
		# acid. Ex. @midi_table['T']['R'] = 62
		def create_midi_table
			# Array of all the amino acids
			@amino_array = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G', '*']

			# Initialize Hashes
			a_hash = Hash.new
			c_hash = Hash.new
			t_hash = Hash.new
			g_hash = Hash.new

			# Fill amino hashes
			@amino_array.each_index{|index|
				a_hash[@amino_array[index]] = decimal_to_hex(index + 44)
				c_hash[@amino_array[index]] = decimal_to_hex(index + 48)
				t_hash[@amino_array[index]] = decimal_to_hex(index + 52)
				g_hash[@amino_array[index]] = decimal_to_hex(index + 56)
			}

			# Place amino hashes in midi_table
			@midi_table['A'] = a_hash
			@midi_table['C'] = c_hash
			@midi_table['T'] = t_hash
			@midi_table['G'] = g_hash
		end

		# Method: dispaly_midi_table
		# => Prints out the midi data in a neat table
		# format. Rows are the nucleic acids; columns are
		# amino acids.
		def display_midi_table
			# write the column names
			@midi_table['A'].each do |key, value|
				print "   #{key}"
			end

			puts # just for a newline

			@midi_table.each do |na, aa_hash|
				print na 
				aa_hash.each do |aa, midi|
					#print "  #{'%02X' % midi}"
					print "  #{midi}"
				end
				puts
			end
		end

		# Method: compute_midi_data
		# => Pulls the midi notes out of midi_table and puts them
		# in an array of arrays.
		def compute_midi_notes

			# Go through each letter of the na sequence and
			# compare it to all possible amino acids.
			# Not currently working because it tries to index
			# nil values. split is returning Fixnums
			index = 0 				# index and aminoIndex make sure that
			aminoIndex = 0 			# the na lines up with the correct aa
			@data.lstrip.split("").each do |i|
				if i =~ /[ACTG]/ # check to see if the char is A C T or G
					if (index % 3) == 0  # which amino the na aligns with depends on which position it is in
				
						puts "#{aminoIndex}  #{@aa_first_pos[aminoIndex]}"
						puts "#{aminoIndex - 1}  #{@aa_second_pos[aminoIndex - 1]}"
						puts "#{aminoIndex - 2}  #{@aa_third_pos[aminoIndex - 2]}"

						if @aa_first_pos[aminoIndex] == nil
							#@midi_notes[0].push(nil)
						else
							@midi_notes[0].push(@midi_table[i][@aa_first_pos[aminoIndex].chr])  
						end

						if @aa_third_pos[aminoIndex - 1] == nil
							#@midi_notes[1].push(nil)
						else
							@midi_notes[1].push(@midi_table[i][@aa_third_pos[aminoIndex - 1].chr])
						end

						if @aa_second_pos[aminoIndex - 1] == nil
							#@midi_notes[2].push(nil)
						else
							@midi_notes[2].push(@midi_table[i][@aa_second_pos[aminoIndex - 1].chr])
						end

					elsif (index % 3) == 1
						puts "#{aminoIndex - 1}  #{@aa_first_pos[aminoIndex - 1]}"
						puts "#{aminoIndex}  #{@aa_second_pos[aminoIndex]}"
						puts "#{aminoIndex + 1}  #{@aa_third_pos[aminoIndex + 1]}"

						if @aa_first_pos[aminoIndex] == nil
							#@midi_notes[0].push(nil)
						else
							@midi_notes[1].push(@midi_table[i][@aa_first_pos[aminoIndex].chr])
						end

						if @aa_second_pos[aminoIndex] == nil
							#@midi_notes[1].push(nil)
						else
							@midi_notes[0].push(@midi_table[i][@aa_second_pos[aminoIndex].chr]) # correct
						end

						if @aa_third_pos[aminoIndex - 1] == nil
							#@midi_notes[2].push(nil)
						else
							@midi_notes[2].push(@midi_table[i][@aa_third_pos[aminoIndex - 1].chr])
						end

					elsif (index % 3) == 2
						puts "#{aminoIndex + 2}  #{@aa_first_pos[aminoIndex + 2]}"
						puts "#{aminoIndex + 1}  #{@aa_second_pos[aminoIndex + 1]}"
						puts "#{aminoIndex}  #{@aa_third_pos[aminoIndex]}"

						if @aa_first_pos[aminoIndex] == nil
							#@midi_notes[0].push(nil)
						else
							@midi_notes[2].push(@midi_table[i][@aa_first_pos[aminoIndex].chr])
						end

						if @aa_second_pos[aminoIndex] == nil
							#@midi_notes[1].push(nil)
						else
							@midi_notes[1].push(@midi_table[i][@aa_second_pos[aminoIndex].chr])
						end

						if @aa_third_pos[aminoIndex] == nil
							#@midi_notes[2].push(nil)
						else
							@midi_notes[0].push(@midi_table[i][@aa_third_pos[aminoIndex].chr]) # correct
						end

						aminoIndex += 1
					end
				end
				index += 1
			end
		end # might need a new algorithm :(

		# Method: display_midi_notes
		# => Prints out the midi notes from the midi_note array.
		# The first row is of notes computed from aa_first_posion
		# and so on. Each column represents the chord derived from
		# the nucleic acids.
		def display_midi_notes
			@midi_notes[0].each{|i| 
				if i == nil 
					print "nil"
				else
					print "#{i} "
				end}
			puts
			@midi_notes[1].each{|i| 
				if i == nil 
					print "nil"
				else
					print "#{i} "
				end}
			puts
			@midi_notes[2].each{|i| 
				if i == nil 
					print "nil"
				else
					print "#{i} "
				end}
			puts # to add a newline
		end

		# Method: create_midi_file(file_name)
		# => Param: file_name - The name of the file to write the 
		# midi data to. Recomended that the file end with the 
		# extension '.mid'
		def create_midi_file(file_name) # needs to be opened with hex encoding
			# Open a file for writing and name it appropriately
			File.open(file_name, 'w') do |midi_file|
				midi_file.print "4D5468640000000600010003#{"%04X" % @tempo}"
				midi_file.print write_track(0)
				midi_file.print write_track(1)
				midi_file.print write_track(2)
			midi_file.close
			end
		end

		# Method: write_midi_header
		# => Writes the header for a standard, type-1 midi file.
		# It does not include the track header. It is overridden
		# in the MidiFormatter class.
		def write_midi_header
			"4D5468640000000600010003#{decimal_to_hex(@tempo)}"
		end # this might be a useless function

		# Method: write_track_header
		# => Writes the header of the result of compute_midi_data.
		def write_track(track)
			this_track = "4D54726B#{calculate_num_bytes(track)}" + @track_data[track]
			this_track
		end # for some reason this is putting a nil at the end

		# Method: create_track_data
		# => Params: track - the index of the track in the midi_notes array
		# must be 0, 1, or 2
		# Return: the data in a track, not including the track header.
		def create_track_data(track)

			@track_data[track] = "" # Initialize it to be a string

			if track == 0
				for i in 0 ... (@midi_notes[track].size - 3)
					@track_data[track] = @track_data[track] + "0090#{@midi_notes[0][i]}60"
					@track_data[track] = @track_data[track] + "#{@tempo}80#{@midi_notes[0][i]}00"
				end
			elsif track == 1
				for i in 1 ... (@midi_notes[track].size - 2)
					@track_data[track] = @track_data[track] + "0090#{@midi_notes[1][i]}60"
					@track_data[track] = @track_data[track] + "#{@tempo}80#{@midi_notes[1][i]}00"
				end
			elsif track == 2 
				for i in 2 ... (@midi_notes[track].size - 1)
					@track_data[track] = @track_data[track] + "0090#{@midi_notes[2][i]}60"
					@track_data[track] = @track_data[track] + "#{@tempo}80#{@midi_notes[2][i]}00"
				end
			end # nil value might be causing holes in the data

			# @midi_notes[track].each do |note|
			# 	# need to figure out a way to start and stop notes
			# 	@track_data[track] = @track_data[track] + "0090#{note}60"
			# 	@track_data[track] = @track_data[track] + "#{@tempo}80#{note}00"
			# end
			@track_data[track] = @track_data[track] + "00FF2F00"
		end 

		# Method: decimal_to_hex
		# => Param: decimal - The decimal number to convert to hex.
		# => Returns @tempo as hexedecimal value properly formatted
		# for the midi header.
		def decimal_to_hex(decimal)
			temp = '%02X' % decimal
			# Make sure an even number of bytes is always written.
			if temp.length.odd?
				'0' + temp
			end
			temp
		end

		# Method: calculate_num_bytes
		# => Param: track - which track to calculate the number of bytes
		# for. Must be a 0, 1, or 2
		# => Return: the number of bytes
		def calculate_num_bytes(track)
			num_bytes = "%08X" % (@track_data[track].length/2) # This is not correct
		end

	end
end

=begin

Notes to self:

Nil values in the note data are causing errors in the midi file.
Easy solution is to not push nil values into midi_notes in the
compute_midi_notes method. This will mean that notes wont align 
as wanted in the output midi file.
Other solution - for every nil value write a dummy note either in the
create_midi_notes method or the create_track_data method. Probably
easier to do in create track data.
Best solution: Find where the nil values are coming from. If everything
is working properly they shouldn't be there.	

Problem might also be in the aminoIndex incrementer.

	
=end