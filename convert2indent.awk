#!/usr/bin/awk -f 

# This script is intended to do the following:
# replace all the fixed width constructs with a 4 level
# indentation. It is automated, so fast, but still requires
# eyeballs to make sure multiline things are readable.
#
# It DOES NOT handle gotos
#
# Furthermore, it ONLY handles line continuation using the '&',
# NOT the line 5 arbitrary character...

# Begin by initializing the indent level to 0
BEGIN{
    il=0;
    isif=0;
    isand=0;
    iselse=0;
}

# If the line begins with an end statement, decrement the indent level
tolower($0)~/^[0-9]*\s*end/{il--;}
tolower($0)~/^[0-9]*\s*else/{il--;iselse=1;}

# reset for subroutines and functions
tolower($0)~/^[0-9]*\s*subroutine/{il=0} 
tolower($0)~/^[0-9]*\s*program/{il=0} 
tolower($0)~/^\s*[a-z]*\s*function/{il=0} 

#Print out the line preceded by its indent count
{
    # remove the leading whitespace
    toprint=$0
    sub(/^\s*/,"",toprint)
    # if not a comment, make it lower case
    if( $0 !~ /^\s*!/ ){
        torpint = tolower(toprint)
    }
    gsub(/\.le\./," <= ",toprint)
    gsub(/\.LE\./," <= ",toprint)
    gsub(/\.lt\./," < ",toprint)
    gsub(/\.LT\./," < ",toprint)
    gsub(/\.ge\./," >= ",toprint)
    gsub(/\.GE\./," >= ",toprint)
    gsub(/\.gt\./," > ",toprint)
    gsub(/\.GT\./," > ",toprint)
    gsub(/\.ne\./," /= ",toprint)
    gsub(/\.NE\./," /= ",toprint)
    gsub(/\.eq\./," == ",toprint)
    gsub(/\.EQ\./," == ",toprint)
    gsub(/\t/,"    ",toprint)
    # prepend correct indent
    for (i=0; i<il; i++){
        toprint="    " toprint
    }
    printf("%s\n",toprint)
    if (iselse){
        il++;
        iselse=0;
    }
}

# Handle multiline statements
/.*&\s*$/ && !isand{
    isand=1;
    il++;
    }
($0 !~ /.*&\s*$/) && isand{
    isand=0;
    il--;
    }

# Handle multiline if statements

# If the line is an if...& statement
tolower($0)~/^[0-9]*\s*if.*&\s*$/{isif=1;}
# if in the middle of an if statement, and ends with a then
tolower($0)~/.*then\s*$/ && isif{
    isif=0;
    il++;
    }
# If in the middle of an if statement, and line doesn't end with a &
/.*[^&]\s*$/ && isif {isif=0;}

# Increment indent level based on subroutine, if...then, do, 
tolower($0)~/^[0-9]*\s*subroutine/{il++} 
tolower($0)~/^[0-9]*\s*program/{il++} 
tolower($0)~/^\s*[a-z]*\s*function/{il++} 
tolower($0)~/^[0-9]*\s*do[^a-z0-9]/{il++} 
tolower($0)~/^[0-9]*\s*if.*then\s*$/{il++;}

# At the end, make sure the indent level is 0
END{
    if(isif != 0 || il != 0){
        print "BAD ENDING";
    }else{
        print "GOOD ENDING";
    }
    }
