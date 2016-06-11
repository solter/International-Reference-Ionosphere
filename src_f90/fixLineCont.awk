#!/usr/bin/awk -f 

# This script will change all line indentation
# from using ampersands at the beginning of a line to the
# end of a previous line

BEGIN{prev_line=0}

{
    if (prev_line){
        # If old style line continuation
        if ( $0 ~ /^\s*&/ && prev_line !~ /.*&\s*$/ ){
            prev_line = prev_line "&"
        }
        printf("%s\n",prev_line)
    }
    sub(/&/," ")
    prev_line = $0
}

END{printf("%s",prev_line)}
