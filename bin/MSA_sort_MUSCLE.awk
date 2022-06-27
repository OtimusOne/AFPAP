/^>/{
    printSeq=0
}

/^>input_sequence/{
    printSeq=1
}

{
    if (printSeq) {
        sub("input_sequence[|]", "", $0)
        print
    } else {
        if (rest) {
            rest = rest "\n" $0
        } else {
            rest = $0
        }
    }
}

END {
    print rest
}
