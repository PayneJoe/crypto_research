#### P2TR

- construct taproot with provided huffman tree leaf scripts

    script consists three parts:
    - covenant script
    - state hash check script
    - state check script

<br />

- construct tweated pubkey with taproot and untweated pubkey $P$

    $$
        Q = P + H(P || c) * G
    $$
    where $Q$ is the tweated pubkey, $c$ is taproot (hash) and $P$ is untweated pubkey.

    we usually use the compressed one, the $x$-coordinate of $Q$ and parity sign

<br />

- construct **P2TR** script pubkey, actually it's just a script, with tweated pubkey $Q$

    ```python3
        [1, Q.x]
    ```
    where `1` is $1$-byte of witness program version for **P2TR**, and $Q.x$ is $32$-bytes of the tweated pubkey

<br />

#### Sequential Covenant

- initial state hash

    <br />

- get script pubkey, tweated pubkey with taproot

    <br />

- construct *caboose* (P2WSH pubkey) of initial transaction

    ```python3
    script = [OP_RETURN, OP_PUSHBYTES_36, init_state_hash(4), init_randomizer(32)]
    script_hash = SHA256(script)
    script_pubkey = [1, script_hash]
    ```
    <br />

- construct init transaction

    ```json
    {
        "input": [
            {
                "previous_output": {
                    "txid": ,
                    "vout": ,
                },
                "script_sig": ,
                "sequence": ,
                "witness": ,
            },
        ],
        "output": [
            {
                "value": 1000000,
                "script_pubkey": ,
            },
            {
                "value": 330,
                "script_pubkey": ,
            }
        ]
    }
    ```
    where `previous_output` specifies the spending input 

    <br />

- compute init transaction id, hash of transaction data

    `witness` is not included

    <br />

- reflection of self-transaction 

    - hash of transaction data $m$

        - epoch, `0`

        - sighash_type with $1$-byte
    
        - tx version with $4$-bytes

        - lock_time with $4$-bytes 

        - spend_type

        - previous_output
            - input $0$
                - prev_txid
                - prev_output index
            - prev_output of input $0$
                - prev_amount
                - script_pubkey locking prev_amount
            - sequence number of input $0$

        - script leaf hash

        - key_version_0, `0`

        - code_separator, `0xFFFFFFFF`

    - generate $e$
        $$
            e = SHA256(TAG || TAG || G || G || m)
        $$

    <br />

    - covenant

        - step $1$, prepare preimage header
            ```python3
            preimage_head = Epoch || HashType || TxVersion || LockTime
            ```
            <br />

        - step $2$, prepare current output1 and dust
            ```python3
            [Pubkey, NewBalance || 34 || PubKey || DUST]
            ```
            <br />

        - step $3$, hash of current outputs
            ```python3
            current_output_hash = SHA256(NewBalance || 34 || PubKey || DUST || 34032 || SHA256(OP_RETURN || OP_PUSHBYTES_36 || NewSateHash || NewRandomizer))

            [Pubkey, OldStateHash, preimage_head || current_output_hash]
            ```
            <br />

        -  step $4$ and $5$, current input1, previous output1, current script hash
            ```python3
            [Pubkey, OldStateHash, PrevBalance, PrevTxid, preimage_head || current_output_hash || SpendType || PrevTxid || 0 || PrevBalance || 34 || Pubkey || Sequence || script_hash || 0 || CodeSeparator]
            ```
            <br />

        - step $6$, compute signature hash $e$ and check with hint
            ```python3
            tag_hash = SHA256(TagHash || TagHash || preimage_head || current_output_hash || SpendType || PrevTxid || 0 || PrevBalance || 34 || Pubkey || Sequence || script_hash || 0 || CodeSeparator)

            sig_hash = SHA256(ChallengeHash || ChallengeHash || G || G || tag_hash)

            [Pubkey, OldStateHash, PrevBalance, PrevTxid, e[:-1], e, sig_hash, OP_EQUALVERIFY]

            [Pubkey, OldStateHash, PrevBalance, PrevTxid, G || e[:-1] || 2129, G, OP_CHECKSIGVERIFY]
            ```
            <br />

- reflection of prev-transaction

    <br />