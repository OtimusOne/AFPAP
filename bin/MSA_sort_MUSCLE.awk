/^>/{
    printSeq=0
}

/^>query_sekvence/{
    printSeq=1
}

{
    if (printSeq) {
        sub("query_sekvence[|]", "", $0)
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
