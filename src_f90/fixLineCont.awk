#!/usr/bin/awk -f 

# This script will change all line indentation
# from using ampersands at the beginning of a line to the
# end of a previous line

BEGIN{prev_line=0}

{
    if (prev_line){
        # If old style line continuation
        if ( $0 ~ /^     [^ ]/ ){
            prev_line = prev_line "&"
            sub(/^     [^ ]/,"      ")
        }else if ( $0 ~ /^\s*&/ && prev_line !~ /.*&\s*$/ ){
            prev_line = prev_line "&"
            sub(/&/," ")
        }
        printf("%s\n",prev_line)
    }
    # strip off trailing characters
    sub(/\s*$/,"",$0)
    prev_line = $0
}

END{printf("%s",prev_line)}
