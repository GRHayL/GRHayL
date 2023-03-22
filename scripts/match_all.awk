# Source: https://unix.stackexchange.com/questions/379385
function match_all(str, regex) {
  out = ""
  start  = RSTART
  end    = RLENGTH
  while( match(str, regex) ) {
    m   = substr(str, RSTART, RLENGTH)
    str = substr(str, RSTART + (RLENGTH ? RLENGTH : 1))
    out = sprintf("%s%s%s ", out, path, m)
    if( str == "" ) break
  }
  RSTART  = start
  RLENGTH = end
  return out
}

set_path == 1 {
  if( curr_file != FILENAME ) {
    curr_file = FILENAME
    path=FILENAME
    sub("make.code.defn", "", path)
  }
}

$0 ~ re_start, $0 ~ re_end {
  matches = sprintf("%s%s", matches, match_all($0, re_in))
}

END {
  print matches
}
