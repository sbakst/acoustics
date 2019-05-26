
input_text_file$ = "/media/sf_2019_reqp/sub_bmps/118/open_me_in_praat.txt"

strings = Read Strings from raw text file: input_text_file$
string_number = Get number of strings

for s to string_number
  selectObject(strings)
  string$ = Get string: s
  Read from file: string$
endfor