#!/usr/bin/awk -f

# Each input line is in the format: sequence_id<space>coverageInPercent<space>identityInPercent
# Output is the sequence_ids that pass the filter.
{
    id = $1
    cover = $2
    ident = $3
    if (cover >= 80 && ident >= 35 && ident <= 95) {
        print id
    }
    next
}
