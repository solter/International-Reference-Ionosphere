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
}

# If the line begins with an end statement, decrement the indent level
tolower($0)~/^[0-9]*\s*end/{il--;}

#Print out the line preceded by its indent count
{printf("%03d+++%s\n",il,$0)}

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
tolower($0)~/^[0-9]*\s*do/{il++} 
tolower($0)~/^[0-9]*\s*if.*then\s*$/{il++;}

# At the end, make sure the indent level is 0
END{
    if(isif != 0 || il != 0){
        print "BAD ENDING";
    }else{
        print "GOOD ENDING";
    }
    }
