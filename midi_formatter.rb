#!/usr/local/bin/ruby
#
# The MidiFormatter class is able to take any number of sequences or fasta
# files and turn them into a single midi file.
#


require 'bio/db/midi'

module Bio
	class MidiFormatter

		attr_reader :midi_format_array

		# Not sure what needs to be done here, yet.
		def initialize

			# The speed of the music
			@tempo = 77

			# The number of MidiFormat objects
			@num_midi_seqs = 0

			# An array of MidiFormat objects
			@midi_format_array = Array.new
		end

		# Method: add_sequence
		# => Param: seq - String na sequence.
		# => Creates a MidiFormat of the sequence and places it
		# in @midi_format_array. For use with sequences not in
		# fasta format.
		def add_sequence(seq)

			# Append false gi
			fasta = ">gi|00000|gb|00000.0|N/A\n" + seq
			
			@midi_format_array.push(Bio::MidiFormat.new(fasta))

			# Increment number of sequences
			@num_midi_seqs += 1
		end # Not currently working

		# Method: add_fasta
		# => Param: fasta - A fasta formatted sequence.
		# => fasta is turned into a MidiFormat objects and stored
		# in @midi_format_array
		def add_fasta(fasta)
			# Create MidiFormat object and push it onto array
			@midi_format_array.push(Bio::MidiFormat.new(fasta))

			# Increment number of sequences
			@num_midi_seqs += 1
		end

		# Method: write_midi_header
		# => 
		def write_midi_header
			"4D546864000000060001#{"%04X" % (@num_midi_seqs*3)}00#{@tempo}"
		end

		def write_midi_data
			@midi_format_array.each do |midi_track|
				midi_track.track_data[0] << midi_track.track_data[1] << midi_track.track_data[2]
			end
		end

		def create_midi_file(file_name)

			# Open the file
			File.open(file_name, 'w') do |midi_file|
				midi_file.print write_midi_header
				@midi_format_array.each do |midi_track|
				midi_file.print midi_track.write_all_tracks
			end
				midi_file.close
			end
		end

	end # end class
end # end module