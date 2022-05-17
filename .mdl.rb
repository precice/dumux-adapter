# Enable all rules by default
all

######################################################
# Deactivated/excluded rules
######################################################


######################################################
# Adapted rules
######################################################

# Enforce dashes as symvol in lists
rule 'MD004', :style => :dash
# Enforce indent with four space characters. This helps with
# mixing ordered and unordered lists
rule 'MD007', :indent => 4
# Extend line length for text.
# This will complain for overly wide tables and code blocks.
rule 'MD013', :line_length => 99999

# Nested lists should be indented with four spaces.
# Modification: Question marks should be allowed
rule 'MD026', :punctuation => '.,;:!'

# Ordered lists must have prefix that increases in order
rule 'MD029', :style => :ordered
