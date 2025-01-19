/*
 * Copyright 2020, Yun (Leo) Zhang <imzhangyun@gmail.com>
 *
 * This file is part of HISAT-3N.
 *
 * HISAT-3N is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT-3N is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT-3N.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef UTILITY_3N_TABLE_H
#define UTILITY_3N_TABLE_H

#include <mutex>
#include <condition_variable>
#include <iostream>
#include <queue>
#include <algorithm>

using namespace std;

/**
 * return complement of input base.
 */
char asc2dnacomp[] = {
        /*   0 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  16 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  32 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,'-',  0,  0,
        /*  48 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  64 */ 0,'T','V','G','H',  0,  0,'C','D',  0,  0,'M',  0,'K','N',  0,
        /*    A   B   C   D           G   H           K       M   N */
        /*  80 */ 0,  0,'Y','S','A',  0,'B','W',  0,'R',  0,  0,  0,  0,  0,  0,
        /*        R   S   T       V   W       Y */
        /*  96 */ 0,'T','V','G','H',  0,  0,'C','D',  0,  0,'M',  0,'K','N',  0,
        /*   a   b   c   d           g   h           k       m   n */
        /* 112 */ 0,  0,'Y','S','A',  0,'B','W',  0,'R',  0,  0,  0,  0,  0,  0,
        /*        r   s   t       v   w       y */
        /* 128 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 144 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 160 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 176 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 192 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 208 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 224 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 240 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

/**
 * the simple data structure to bind quality score and position (on reference) together.
 */
class PosQuality {
public:
    int readPos; // 0-based
    int refPos;  // 0-based
    char qual;
    bool converted;
    bool remove;

    PosQuality(int& inputPos) {
        readPos = inputPos;
        refPos = inputPos;
        remove = true;
    }

    void setQual (char& inputQual, bool inputConverted) {
        qual = inputQual;
        converted = inputConverted;
        remove = false;
    }
};

/**
 * the base class for string we need to search.
 */
class string_search {
public:
    int start;
    string s;
    int stringLen;

    void initialize() {
        start = 0;
        stringLen = 0;
        s.clear();
    }

    void loadString(string intputString) {
        s = intputString;
        stringLen = s.size();
        start = 0;
    }
};


/**
 * to store CIGAR string and search segments in it.
 */
class CIGAR : public string_search{
public:

    bool getNextSegment(int& len, char& symbol) {
        if (start == stringLen) {
            return false;
        }
        len = 0;
        int currentIndex = start;
        while (true) {
            if (isalpha(s[currentIndex])) {
                len = stoi(s.substr(start, currentIndex-start));
                symbol = s[currentIndex];
                start = currentIndex+1;
                return true;
            }
            currentIndex++;
        }
    }
};

/**
 * to store MD tag and search segments in it.
 */
class MD_tag : public string_search {
public:

    bool getNextSegment(string& seg) {
        if (start >= stringLen) {
            return false;
        }
        seg.clear();
        int currentIndex = start;
        bool deletion = false;

        while (true) {
            if (currentIndex >= stringLen) {
                start = currentIndex + 1;
                return !seg.empty();
            }
            if (seg.empty() && s[currentIndex] == '0') {
                currentIndex++;
                continue;
            }
            if (isalpha(s[currentIndex])) {
                if (seg.empty()) {
                    seg = s[currentIndex];
                    start = currentIndex+1;
                    return true;
                } else {
                    if (deletion) {
                        seg += s[currentIndex];
                        //currentIndex++;
                    } else {
                        start = currentIndex;
                        return true;
                    }
                }
            } else if (s[currentIndex] == '^') {
                if (seg.empty()) {
                    seg = s[currentIndex];
                    deletion = true;
                } else {
                    start = currentIndex;
                    return true;
                }
            } else { // number
                if (seg.empty()) {
                    seg = s[currentIndex];
                } else {
                    if (deletion || isalpha(seg.back())) {
                        start = currentIndex;
                        return true;
                    } else {
                        seg += s[currentIndex];
                    }
                }
            }
            currentIndex++;
        }
    }
};

/**
 * simple safe queue
 */
template <typename T>
class SafeQueue {
private:
    mutex mutex_;
    queue<T> queue_;
    std::condition_variable notEmpty_;

    string getReadName(string* line){
        int startPosition = 0;
        int endPosition;

        endPosition = line->find("\t", startPosition);
        string readName = line->substr(startPosition, endPosition - startPosition);
        return readName;
    }

public:
    void pop() {
        std::unique_lock<std::mutex> lk(mutex_);
        queue_.pop();
    }

    T front() {
        std::unique_lock<std::mutex> lk(mutex_);
        T value = queue_.front();
        return value;
    }

    int size() {
        std::unique_lock<std::mutex> lk(mutex_);
        int s = queue_.size();
        return s;
    }

    /**
     * return true if the queue is not empty and pop front and get value.
     * return false if the queue is empty.
     */
    bool popFront(T& value) {
        std::unique_lock<std::mutex> lk(mutex_);
        bool isEmpty = queue_.empty();
        if (!isEmpty) {
            value = queue_.front();
            queue_.pop();
        }
        return !isEmpty;
    }

    void push(T value) {
        std::unique_lock<std::mutex> lk(mutex_);
        queue_.push(value);
    }

    bool empty() {
        std::unique_lock<std::mutex> lk(mutex_);
        bool check = queue_.empty();
        return check;
    }

    void pushAndNotify(T value) {
try {
        // outputPositionPool get something to run, thus notify other to wake
        std::unique_lock<std::mutex> lk(mutex_);
        if (queue_.empty()) {
            queue_.push(value);
            std::cerr << "HEY. YOU SHOULD WAKE UP\t";
            notEmpty_.notify_all();
        } else {
            queue_.push(value);
        }
        lk.unlock();
} catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << std::endl;
    exit(1);
}
    }

    // You should have done it as a virtual base(i.e. interface), but I have no time reconstruct it
    bool printOrWait(T& pos) {
        std::unique_lock<std::mutex> lk(mutex_);
        if (queue_.size()) {
            pos = queue_.front();
            lk.unlock();
            return true;
        } else {
            // this_thread::sleep_for (std::chrono::microseconds(1));
            std::cerr << "I wait for work??\t";
            notEmpty_.wait(lk);
            lk.unlock();
            std::cerr << "I get waken up\t";
            return false;
        }
    }
};

/**
 * store one chromosome and it's stream position
 */
class ChromosomeFilePosition {
public:
    string chromosome;
    streampos linePos;
    ChromosomeFilePosition(string inputChromosome, streampos inputPos) {
        chromosome = inputChromosome;
        linePos = inputPos;
    }

    bool operator < (const ChromosomeFilePosition& in) const{
        return chromosome < in.chromosome;
    }
};

/**
 * store all chromosome and it's stream position
 */
class ChromosomeFilePositions {
public:
    vector <ChromosomeFilePosition> pos;

    /**
     * input the chromosome name and it's streamPos, if it is not in pos, add it.
     */
    void append (string &chromosome, streampos& linePos) {
        pos.push_back(ChromosomeFilePosition(chromosome, linePos));
    }

    /**
     * make binary search on pos for target chromosome name
     */
    int findChromosome(string &targetChromosome, int start, int end) {
        if (start <= end) {
            int middle = (start + end) / 2;
            if (pos[middle].chromosome == targetChromosome) {
                return middle;
            }
            if (pos[middle].chromosome > targetChromosome) {
                return findChromosome(targetChromosome, start, middle-1);
            }
            return findChromosome(targetChromosome, middle+1, end);
        }
        else
        {
            // cannot find the chromosome! throw!
            cerr << "Cannot find the chromosome: " << targetChromosome << " in reference file." << endl;
            throw 1;
        }
    }

    /**
     * given targetChromosome name, return its streampos
     */
    streampos getChromosomePosInRefFile(string &targetChromosome)
    {
        int index = findChromosome(targetChromosome, 0, pos.size()-1);
        assert(pos[index].chromosome == targetChromosome);
        return pos[index].linePos;
    }

    /**
     * sort the pos by chromosome name
     */
    void sort()
    {
        std::sort(pos.begin(), pos.end());
    }
};
#endif //UTILITY_3N_TABLE_H
