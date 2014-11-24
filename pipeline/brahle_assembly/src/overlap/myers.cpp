#include "myers.h"

#include <stdint.h>

namespace overlap {

namespace {

typedef uint64_t Word;
const int32_t kWordSize = sizeof(Word) * 8;
const Word kHighBitMask = ((Word)1) << (kWordSize-1);

/**
 * Corresponds to Advance_Block function from Myers.
 * Calculates one word(block), which is part of a column.
 * Highest bit of word is most bottom cell of block from column.
 * @param [in] Pv  Bitset, Pv[i] == 1 if vin is +1, otherwise Pv[i] == 0.
 * @param [in] Mv  Bitset, Mv[i] == 1 if vin is -1, otherwise Mv[i] == 0.
 * @param [in] Eq  Bitset, Eq[i] == 1 if match, 0 if mismatch.
 * @param [in] hin  Will be +1, 0 or -1.
 * @param [out] PvOut  Bitset, PvOut[i] == 1 if vout is +1, otherwise PvOut[i] == 0.
 * @param [out] MvOut  Bitset, MvOut[i] == 1 if vout is -1, otherwise MvOut[i] == 0.
 * @param [out] hout  Will be +1, 0 or -1.
 */
inline int calculateBlock(Word Pv, Word Mv, Word Eq, const int hin,
                                 Word &PvOut, Word &MvOut) {
    Word Xv = Eq | Mv;
    if (hin < 0)
        Eq |= (Word)1;
    Word Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;

    Word Ph = Mv | ~(Xh | Pv);
    Word Mh = Pv & Xh;

    int hout = 0;
    if (Ph & kHighBitMask)
        hout = 1;
    else if (Mh & kHighBitMask)
        hout = -1;

    Ph <<= 1;
    Mh <<= 1;

    if (hin < 0)
        Mh |= (Word)1;
    else if (hin > 0)
        Ph |= (Word)1;
    PvOut = Mh | ~(Xv | Ph);
    MvOut = Ph & Xv;

    return hout;
}

/**
 * Does ceiling division x / y.
 * Note: x and y must be non-negative and x + y must not overflow.
 */
inline int ceilDiv(int x, int y) {
    return (x + y - 1) / y;
}

inline int min(int x, int y) {
    return x < y ? x : y;
}

inline int max(int x, int y) {
    return x > y ? x : y;
}

inline int abs(int x) {
    return x < 0 ? -1 * x : x;
}

/**
 * @param [in] mode  MYERS_MODE_HW or MYERS_MODE_SHW
 */
int myersCalcEditDistanceSemiGlobal(Word* P, Word* M, int* score, Word** Peq, int W, int maxNumBlocks,
                                           const unsigned char* target, int targetLength,
                                           int alphabetLength, int k, int mode, int* bestScore_, int* position_) {
    // firstBlock is 0-based index of first block in Ukkonen band.
    // lastBlock is 0-based index of block AFTER last block in Ukkonen band. <- WATCH OUT!
    int firstBlock = 0;
    int lastBlock = min(ceilDiv(k + 1, kWordSize), maxNumBlocks); // y in Myers

    // Initialize P, M and score
    for (int b = 0; b < lastBlock; b++) {
        score[b] = (b+1) * kWordSize;
        P[b] = (Word)-1; // All 1s
        M[b] = (Word)0;
    }

    int bestScore = -1;
    int position = -1;
    for (int c = 0; c < targetLength + W; c++) { // for each column
        // We pretend like target is padded at end with W wildcard symbols
        Word* Peq_c = c < targetLength ? Peq[target[c]] : Peq[alphabetLength];

        //----------------------- Calculate column -------------------------//
        int hout = mode == MYERS_MODE_HW ? 0 : 1;
        for (int b = firstBlock; b < lastBlock; b++) {
            hout = calculateBlock(P[b], M[b], Peq_c[b], hout, P[b], M[b]);
            score[b] += hout;
        }
        //------------------------------------------------------------------//

        //---------- Adjust number of blocks according to Ukkonen ----------//
        if (score[firstBlock] >= k + kWordSize)
            firstBlock++;

        if ((score[lastBlock-1] - hout <= k) && (lastBlock < maxNumBlocks)
            && ((Peq_c[lastBlock] & (Word)1) || hout < 0)) {
            // If score of left block is not too big, calculate one more block
            lastBlock++;
            int b = lastBlock-1; // index of last block (one we just added)
            P[b] = (Word)-1; // All 1s
            M[b] = (Word)0;
            score[b] = score[b-1] - hout + kWordSize + calculateBlock(P[b], M[b], Peq_c[b], hout, P[b], M[b]);
        }
        else
            while (lastBlock > 0 && score[lastBlock-1] >= k + kWordSize)
                lastBlock--;

        // If band stops to exist, return -1
        if (lastBlock <= firstBlock) {
            *bestScore_ = -1;
            *position_ = -1;
            return MYERS_STATUS_OK;
        }
        //------------------------------------------------------------------//

        //------------------------- Update best score ----------------------//
        if (c >= W && lastBlock == maxNumBlocks) { // We ignore scores from first W columns, they are not relevant.
            int colScore = score[maxNumBlocks-1];
            if (colScore <= k) { // Scores > k dont have correct values (so we cannot use them), but are certainly > k.
                // NOTE: Score that I find in column c is actually score from column c-W
                if (bestScore == -1 || colScore < bestScore) {
                    bestScore = colScore;
                    position = c - W;
                }
            }
        }
        //------------------------------------------------------------------//
    }

    *bestScore_ = bestScore;
    *position_ = position;
    return MYERS_STATUS_OK;
}





int myersCalcEditDistanceNW(Word* P, Word* M, int* score, Word** Peq, int W, int maxNumBlocks,
                                   int queryLength,
                                   const unsigned char* target, int targetLength,
                                   int k, int* bestScore_, int* position_) {

    if (k < abs(targetLength - queryLength)) {
        *bestScore_ = *position_ = -1;
        return MYERS_STATUS_OK;
    }

    // firstBlock is 0-based index of first block in Ukkonen band.
    // lastBlock is 0-based index of block AFTER last block in Ukkonen band. <- WATCH OUT!
    int firstBlock = 0;
    // Added + 1 below, without it about 2 of 100000 test examples fail
    int lastBlock = min(ceilDiv(k - abs(targetLength - queryLength) + 1, kWordSize), maxNumBlocks); // y in Myers

    // Initialize P, M and score
    for (int b = 0; b < lastBlock; b++) {
        score[b] = (b+1) * kWordSize;
        P[b] = (Word)-1; // All 1s
        M[b] = (Word)0;
    }

    for (int c = 0; c < targetLength; c++) { // for each column
        Word* Peq_c = Peq[target[c]];

        //----------------------- Calculate column -------------------------//
        int hout = 1;
        for (int b = firstBlock; b < lastBlock; b++) {
            hout = calculateBlock(P[b], M[b], Peq_c[b], hout, P[b], M[b]);
            score[b] += hout;
        }
        //------------------------------------------------------------------//

        //---------- Adjust number of blocks according to Ukkonen ----------//
        // Adjust first block
        if (score[firstBlock] >= k + kWordSize) // TODO: put some stronger constraint
            firstBlock++;

        // Adjust last block
        if (score[lastBlock-1] - hout // score of block to left
            + max(0, targetLength - c/*column of block to left*/)
            + max(0, queryLength - (lastBlock * kWordSize) <= k)
            && (lastBlock < maxNumBlocks)
            && ((Peq_c[lastBlock] & (Word)1) || hout < 0)) {
            // If score of left block is not too big, calculate one more block
            lastBlock++;
            int b = lastBlock-1; // index of last block (one we just added)
            P[b] = (Word)-1; // All 1s
            M[b] = (Word)0;
            score[b] = score[b-1] - hout + kWordSize + calculateBlock(P[b], M[b], Peq_c[b], hout, P[b], M[b]);
        } else {
            while (lastBlock > 0
                   && score[lastBlock-1] >= k - max(0, targetLength - (c + 1))
                   - max(0, queryLength - (lastBlock * kWordSize)) + kWordSize) {
                lastBlock--;
            }
        }

        // If band stops to exist, return -1
        if (lastBlock <= firstBlock) {
            *bestScore_ = *position_ = -1;
            return MYERS_STATUS_OK;
        }
        //------------------------------------------------------------------//

    }

    if (lastBlock == maxNumBlocks) { // If last block of last column was calculated
        int bestScore = score[maxNumBlocks-1];

        for (int i = 0; i < W; i++) {
            if (P[maxNumBlocks-1] & kHighBitMask)
                bestScore--;
            if (M[maxNumBlocks-1] & kHighBitMask)
                bestScore++;
            P[maxNumBlocks-1] <<= 1;
            M[maxNumBlocks-1] <<= 1;
        }

        if (bestScore <= k) {
            *bestScore_ = bestScore;
            *position_ = targetLength - 1;
            return MYERS_STATUS_OK;
        }
    }


    *bestScore_ = *position_ = -1;
    return MYERS_STATUS_OK;
}

}  // unnamed namespace

int MyersEditDistance(const unsigned char* query, int queryLength,
                      const unsigned char* target, int targetLength,
                      int alphabetLength, int k, int mode, int* bestScore, int* position) {

    /*--------------------- INITIALIZATION ------------------*/
    int maxNumBlocks = ceilDiv(queryLength, kWordSize); // bmax in Myers

    Word* P = new Word[maxNumBlocks]; // Contains Pvin for each block (column is divided into blocks)
    Word* M = new Word[maxNumBlocks]; // Contains Mvin for each block
    int* score = new int[maxNumBlocks]; // Contains score for each block
    Word** Peq = new Word*[alphabetLength+1]; // [alphabetLength+1][maxNumBlocks]. Last symbol is wildcard.

    int W = maxNumBlocks * kWordSize - queryLength; // number of redundant cells in last level blocks

    // Build Peq (1 is match, 0 is mismatch). NOTE: last column is wildcard(symbol that matches anything) with just 1s
    for (int symbol = 0; symbol <= alphabetLength; symbol++) {
        Peq[symbol] = new Word[maxNumBlocks];
        for (int b = 0; b < maxNumBlocks; b++) {
            if (symbol < alphabetLength) {
                Peq[symbol][b] = 0;
                for (int r = (b+1) * kWordSize - 1; r >= b * kWordSize; r--) {
                    Peq[symbol][b] <<= 1;
                    // NOTE: We pretend like query is padded at the end with W wildcard symbols
                    if (r >= queryLength || query[r] == symbol)
                        Peq[symbol][b] += 1;
                }
            } else { // Last symbol is wildcard, so it is all 1s
                Peq[symbol][b] = (Word)-1;
            }
        }
    }
    /*-------------------------------------------------------*/


    /*------------------ MAIN CALCULATION -------------------*/
    *bestScore = -1;
    *position = -1;
    if (k < 0) { // If valid k is not given, auto-adjust k until solution is found.
        k = kWordSize;
        while (*bestScore == -1) {
            if (mode == MYERS_MODE_HW || mode == MYERS_MODE_SHW)
                myersCalcEditDistanceSemiGlobal(P, M, score, Peq, W, maxNumBlocks,
                                                target, targetLength,
                                                alphabetLength, k, mode, bestScore, position);
            else  // mode == MYERS_MODE_NW
                myersCalcEditDistanceNW(P, M, score, Peq, W, maxNumBlocks,
                                        queryLength, target, targetLength,
                                        k, bestScore, position);
            k *= 2;
        }
    } else {
        if (mode == MYERS_MODE_HW || mode == MYERS_MODE_SHW)
            myersCalcEditDistanceSemiGlobal(P, M, score, Peq, W, maxNumBlocks,
                                            target, targetLength,
                                            alphabetLength, k, mode, bestScore, position);
        else  // mode == MYERS_MODE_NW
            myersCalcEditDistanceNW(P, M, score, Peq, W, maxNumBlocks,
                                    queryLength, target, targetLength,
                                    k, bestScore, position);
    }
    /*-------------------------------------------------------*/

    //--- Free memory ---//
    delete[] P;
    delete[] M;
    delete[] score;
    for (int i = 0; i < alphabetLength+1; i++)
        delete[] Peq[i];
    delete[] Peq;
    //-------------------//

    return MYERS_STATUS_OK;
}


}  // namespace overlap
