#!/usr/bin/awk -f

# Each line is in this format: id<space>coverageInPercent<space>identityInPercent
# Output is only sequence_id that pass the filter.
{
    id = $1
    cover = $2
    ident = $3
    if (cover >= 80 && ident >= 35 && ident <= 95) {
        print id
    }
    next
}
