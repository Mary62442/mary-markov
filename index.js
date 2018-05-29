const findSequence = (sequence, states) => {
    return sequence
    .reduce((all, curr) => {        
        all.push(states.findIndex(x => x.state === curr));        
        return all;        
    }, [])
};

const findIndex = (arr, el) => {
    return arr.indexOf(el);
};

const calculateProb = (stateTrans2, init, sequenceArr) => ({
    sequenceProb : () => {
        return sequenceArr
        .reduce((total, curr, i, arr) => {   
        if (i === 0) total += init[curr];    
        else  total *= stateTrans2[arr[i-1]][curr];         
        return total;    
        },0);
    }    
});

const Bayes = (hmm) => ({
    bayesTheorem : (ob, hState) => {
        let hStateIndex = findIndex(hmm.states, hState);
        let obIndex = findIndex(hmm.observables, ob);
        let emissionProb = hmm.emissionMatrix[obIndex][hStateIndex];
        let initHState = hmm.initialProb[hStateIndex];
        let obProb = hmm.emissionMatrix[obIndex].reduce((total, em, i) => {
            total += (em*hmm.initialProb[i]);
            return total;
        }, 0);

        let bayesResult = (emissionProb*initHState)/obProb;
        return bayesResult;
    }
});

const Viterbi = (hmm) => ({

    viterbiAlgorithm : function(obSequence) {

        //Find first trellis
        let trellis1 = this.firstTrellis(obSequence);

        // Find first hidden state given trellis1
        let firstState = {prob: Math.max(...trellis1), index: findIndex(trellis1, Math.max(...trellis1))};

        //Find the whole hidden state sequence by calling the nextTrellis function.
        let stateSequenceObj = this.nextTrellis(obSequence, trellis1, 1, [firstState]);

        // Find the total probability of the sequence
        let highestProb = stateSequenceObj.reduce((total, curr, i, arr) => {
            if ( i === 0) total+=curr.prob;
            else total*=curr.prob;            
            return total;
        }, 0);

        // Return an object with the state sequence and the total probability
        let stateSequence = stateSequenceObj.map((s,i) => hmm.states[s.index] )
        return {states: stateSequence, prob:highestProb, statesTrellis: stateSequenceObj.map(s => s.prob)};
    },

    // The first trellis is found by multiplying the initial probability of a given
    // state by its emission probability of the first observable in the sequence   

    firstTrellis : function(obSequence) {
        let trellis1 = [];

        // Here we find the index of the first observable in the sequence
        let obIndex = findIndex(hmm.observables, obSequence[0]);

        // Here we find the array of emission probabilities of the observable
        let obEmission = hmm.emissionMatrix[obIndex];       

        // For each emission we multiply its state's corresponding initial probability
        hmm.initialProb.forEach((p,i) => {
            trellis1.push(p*obEmission[i]);
        });

        // Trellis1 represents the state trellises of the first step of the sequence.
        return trellis1;
    },

    // To apply the Viterbi algorithm to the following steps of the sequence we need to
    // provide the observation sequence, the previous trellis array, the step in the observation
    // sequence and the array of previous hidden states already found.

    nextTrellis : function(obSequence, prevTrellis, i, sSequence) {   
        let obIndex = i;

        // If obIndex is equal to the observation sequence length, then we can
        // return the hidden state sequence that we have found.
        if (obIndex === obSequence.length) return sSequence;

        // Else we create for each state trellises of the step obIndex an array which will
        // contain all the Viterbi calculations derived from the previous step's trellises.        
        let nextTrellis = [];
        for (let s = 0; s < hmm.states.length; s++) {
            let trellisArr = [];
            prevTrellis.forEach((prob, i) => {
                let trans = hmm.transMatrix[i][s];
                let emiss = hmm.emissionMatrix[findIndex(hmm.observables, obSequence[obIndex])][s];

                // Viterbi algorithm: previous Viterbi path * transition probability of previous state to
                // current state * emission probability of current observation given current state
                trellisArr.push(prob*trans*emiss);
            });  
            
            // To find the trellis for the next step in the sequence we look at the maximum
            // probability in each of the trellis arrays
            nextTrellis.push(Math.max(...trellisArr));
        };

        // Before calling the nextTrellis function recursively we need to add to the hidden state sequence the
        // highest element of the trellis we have found and record its index so as to find which state it refers to.
        let highestTrellis = {prob: Math.max(...nextTrellis), index: findIndex(nextTrellis, Math.max(...nextTrellis)) };
        sSequence.push(highestTrellis);

        // As long as there are observations in the sequence we call nextTrellis recursively passing the observation sequence as
        // always, the new trellis we have found, the next observation index, and the enriched hidden state sequence.
        return this.nextTrellis(obSequence, nextTrellis, obIndex+1, sSequence);
    }
});

const MarkovChain = (states, init, sequence) => {
    let info = {  
        states: states.map(s => s.state),    
        transMatrix : states.map(s => s.prob),   
        initialProb : init,         
        sequence,
        sequenceArr : findSequence(sequence, states) 
    }     
    return Object.assign({}, info, calculateProb(info.transMatrix, init, info.sequenceArr))
};

const HMM = (states, observables, init) => {
    let hmm = {  
        states: states.map(s => s.state),        
        transMatrix : states.map(s => s.prob),  
        initialProb : init,   
        observables : observables.map( o => o.obs ),
        emissionMatrix : observables.map(o => o.prob)     
    }     
    return Object.assign({}, hmm, Bayes(hmm), Viterbi(hmm))
};

exports.MarkovChain = MarkovChain;
exports.HMM = HMM;

