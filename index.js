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

const gamma = (alpha,beta,forward) => {
    return (alpha*beta)/forward;
};

const xi = (alpha, trans, emiss, beta, forward) => {
    return (alpha*trans*emiss*beta)/forward;
}

const calculateProb = (stateTrans2, init, states) => ({
    sequenceProb : (sequence) => {
        return findSequence(sequence, states) 
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

const Forward = (hmm) => ({
    forwardAlgorithm : function(obSequence) {
        let initAlphas = this.initForward(obSequence);
        let allAlphas = this.recForward(obSequence, initAlphas, 1, [initAlphas]);
        return {alphas: allAlphas, alphaF : this.termForward(allAlphas)};
    },

    initForward : function(obSequence) {
        let initTrellis = [];
        let obIndex = findIndex(hmm.observables, obSequence[0]);
        let obEmission = hmm.emissionMatrix[obIndex];  
        hmm.initialProb.forEach((p,i) => {
            initTrellis.push(p*obEmission[i]);
        });
        return initTrellis;
    },

    recForward : function(obSequence, prevTrellis, i, alphas) {   
        let obIndex = i;        
        if (obIndex === obSequence.length) return alphas;             
        let nextTrellis = [];
        for (let s = 0; s < hmm.states.length; s++) {
            let trellisArr = [];
            prevTrellis.forEach((prob, i) => {
                let trans = hmm.transMatrix[i][s];
                let emiss = hmm.emissionMatrix[findIndex(hmm.observables, obSequence[obIndex])][s];                
                trellisArr.push(prob*trans*emiss);
            });              
            nextTrellis.push(trellisArr.reduce((tot,curr) => tot+curr));
        };      
        alphas.push(nextTrellis);
        return this.recForward(obSequence, nextTrellis, obIndex+1, alphas);
    },

    termForward : function(alphas) {
        return alphas[alphas.length-1]        
        .reduce((tot,val) => tot+val);
    }
});

const Backward = (hmm) => ({
    backwardAlgorithm : function(obSequence) {
        let initBetas = hmm.states.map(s => 1);
        let allBetas = this.recBackward(obSequence, initBetas, obSequence.length-1, [initBetas]);
        return {betas: allBetas, betaF:this.termBackward(obSequence, allBetas)};
    },   
    recBackward : function(obSequence, prevBetas, i, betas) {   
        let obIndex = i;        
        if (obIndex === 0) return betas;             
        let nextTrellis = [];
        for (let s = 0; s < hmm.states.length; s++) {
            let trellisArr = [];
            prevBetas.forEach((prob, i) => {
                let trans = hmm.transMatrix[s][i];
                let emiss = hmm.emissionMatrix[findIndex(hmm.observables, obSequence[obIndex])][i];                
                trellisArr.push(prob*trans*emiss);
            });              
            nextTrellis.push(trellisArr.reduce((tot,curr) => tot+curr));
        };      
        betas.push(nextTrellis);
        return this.recBackward(obSequence, nextTrellis, obIndex-1, betas);
    },

    termBackward : function(obSequence, betas) {
        let finalBetas = betas[betas.length-1].reduce((tot,curr,i) => {
            let obIndex = findIndex(hmm.observables, obSequence[0]);
            let obEmission = hmm.emissionMatrix[obIndex];  
            tot.push(curr*hmm.initialProb[i]*obEmission[i]);
            return tot;
        },[]);
        return finalBetas.reduce((tot,val) => tot+val);
    }
});

const EM = (hmm, forwardObj, backwardBetas, obSequence) => ({
    initialGamma : (stateI) => {
        return gamma(forwardObj.alphas[0][stateI], backwardBetas[0][stateI], forwardObj.alphaF);
    },

    gammaTimesInState : (stateI) => {        
        let gammas = [];
        for ( let t = 0; t < obSequence.length; t++) {
            gammas.push(gamma(forwardObj.alphas[t][stateI], backwardBetas[t][stateI], forwardObj.alphaF));
        };
        return gammas.reduce((tot,curr) => tot+curr);
    },

    gammaTransFromState : (stateI) => {        
        let gammas = [];
        for ( let t = 0; t < (obSequence.length-1); t++) {
            gammas.push(gamma(forwardObj.alphas[t][stateI], backwardBetas[t][stateI], forwardObj.alphaF));
        };
        return gammas.reduce((tot,curr) => tot+curr);
    },

    xiTransFromTo : (stateI, stateJ) => {       
        let xis = [];
        for ( let t = 0; t < (obSequence.length-1); t++) {
            let alpha = forwardObj.alphas[t][stateI];
            let trans = hmm.transMatrix[stateI][stateJ];
            let emiss = hmm.emissionMatrix[findIndex(hmm.observables, obSequence[t+1])][stateJ];            
            let beta = backwardBetas[t+1][stateJ];
            xis.push(xi(alpha,trans,emiss,beta,forwardObj.alphaF));
        };
        return xis.reduce((tot,curr) => tot+curr);
    },

    gammaTimesInStateWithOb : (stateI, obIndex) => {
        let obsK = hmm.observables[obIndex];
        let stepsWithOb = obSequence.reduce((tot,curr,i)=> {
            if (curr === obsK) tot.push(i);
            return tot;
        },[]);
        
        let gammas = [];
        stepsWithOb.forEach( step => {
            gammas.push(gamma(forwardObj.alphas[step][stateI], backwardBetas[step][stateI], forwardObj.alphaF));
        });    
        return gammas.reduce((tot,curr) => tot+curr);
    }
});

const BaumWelch = (hmm) => ({

    baumWelchAlgorithm : (obSequence) => {
        let forwardObj = Forward(hmm).forwardAlgorithm(obSequence);
        let backwardBetas = Backward(hmm).backwardAlgorithm(obSequence).betas.reverse();

        let EMSteps = EM(hmm, forwardObj, backwardBetas, obSequence);

        let initProb = [];
        let transMatrix = [];
        let emissMatrix = [];

        for (let i = 0; i< hmm.states.length; i++) {
            initProb.push(EMSteps.initialGamma(i));
            let stateTrans = [];
            for (let j = 0; j< hmm.states.length; j++) {
                stateTrans.push(EMSteps.xiTransFromTo(i,j)/ EMSteps.gammaTransFromState(i));
            };
            transMatrix.push(stateTrans);
        };

        for (let o = 0; o < hmm.observables.length; o++) {
            let obsEmiss = [];
            for (let i = 0; i< hmm.states.length; i++) {
                obsEmiss.push(EMSteps.gammaTimesInStateWithOb(i,o)/ EMSteps.gammaTimesInState(i));
            };
            emissMatrix.push(obsEmiss);
        };

        let hiddenStates = transMatrix
        .reduce((tot,curr,i) => {
            let stateObj = {state: hmm.states[i], prob: curr}
            tot.push(stateObj);
            return tot;
        }, []);

        let observables = emissMatrix
        .reduce((tot,curr,i) => {
            let obsObj = {obs:hmm.observables[i], prob:curr};
            tot.push(obsObj);
            return tot;
        }, []);

        return HMM(hiddenStates, observables, initProb);
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

const MarkovChain = (states, init) => {
    let info = {  
        states: states.map(s => s.state),    
        transMatrix : states.map(s => s.prob),   
        initialProb : init  
    }     
    return Object.assign({}, info, calculateProb(info.transMatrix, init, states))
};

const HMM = (states, observables, init) => {
    let hmm = {  
        states: states.map(s => s.state),        
        transMatrix : states.map(s => s.prob),  
        initialProb : init,   
        observables : observables.map( o => o.obs ),
        emissionMatrix : observables.map(o => o.prob)     
    }     
    return Object.assign({}, hmm, Bayes(hmm), Viterbi(hmm), Forward(hmm), Backward(hmm), BaumWelch(hmm))
};

exports.MarkovChain = MarkovChain;
exports.HMM = HMM;

