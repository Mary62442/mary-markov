# mary-markov

## An npm package to calculate probabilities from Markov Chains and Hidden Markov Models
  
`const markov = require('mary-markov');`  

## Markov Chain

A simple Markov Chain requires states, transition probabilities, initial probabilities and a state sequence.

The states array must be an array of objects with the properties state and prob.

The prob property is an array representing the corresponding line of the matrix of the state.

For example, given a series of states *S = { 'sunny', 'cloudy', 'rainy'}* the transition matrix would be: 

		| 0.4 0.4 0.2 |	
	A = | 0.3 0.3 0.4 |
		| 0.2 0.5 0.3 |
(represents the transition probabilities between the weather sunny, cloudy and rainy)

We'd instantiate a series of states as such:  

    let states = [
	    {state: 'sunny', prob:[0.4, 0.4, 0.2]},
	    {state: 'cloudy', prob:[0.3, 0.3, 0.4]},
	    {state: 'rainy', prob:[0.2, 0.5, 0.3]}
    ];

The initial probabilities of the states will be a simple array.
Each element of the array has the same index of the corresponding state in the states array.

Therefore, in this example, 0.4 is the sunny probability, 0.3 is the cloudy probability, and the final 0.3 is the rainy probability.

`let init = [0.4, 0.3, 0.3];`

The state sequence is an array of states which will be used to calculate sequence probability.  

`let stateSeq = ['sunny', 'rainy', 'sunny', 'sunny', 'cloudy'];`

To instantiate the Markov Chain we pass the states, the initial probabilities and the state sequence as parameters of the MarkovChain function.

`let markovChain = markov.MarkovChain(states, init, stateSeq);`

### Markov Chain sequence probability
To then calculate the sequence probability we call the sequenceProb() function on the object just instantiated.

	let seqProbability = markovChain.sequenceProb(); //0.002560000000000001

### Markov Chain properties
|Property | Description|
|------------ | -------------|
|states | Array of the names of the states|
|transMatrix | Array of arrays representing thetransition probabilities|
|initialProb | Array of initial probabilities|
|sequence | Sequence array of the input|
|sequenceArr | Index representation of the sequence array compared to the states array|

Example:

    console.log(markovChain.transMatrix) // [ [ 0.4, 0.4, 0.2 ], [ 0.3, 0.3, 0.4 ], [ 0.2, 0.5, 0.3 ] ]

   

## Hidden Markov Model

A Hidden Markov Model requires hidden states, transition probabilities, observables, emission probabilities, and initial probabilities.

For example, given a series of states *S = { 'AT-rich', 'CG-rich'}* the transition matrix would look like this:

		| 0.95 0.05 |
	A = | 			|
		| 0.1  0.9	|
(represents the transition probabilities between AT-rich and CG-rich segments in a DNA sequence)

In the program we'd instantiate a series of hidden states as such:

    let hiddenStates = [    
	    {state: 'AT-rich', prob: [0.95, 0.05]},    
	    {state: 'CG-rich', prob: [0.1, 0.9]}     
    ];

The hidden states array must be an array of objects with the properties state and prob.

The prob property is the array representing the corresponding line of the matrix of the hidden state.


The observables array is similar to the hiddenStates array.
Given a series of observables *O = { 'A', 'T', 'C', 'G' }* the emission probabilities would be represented in the matrix:

		| 0.4  0.1 0.1 0.4 	|
	B = | 					|
		| 0.05 0.4 0.4 0.05 |

(represents the emission probabilities of the observables A, T, C, G given the hidden states AT-rich and CG-rich)

In the program the observables would be instantiated as:

    let observables = [    
	    {obs: 'A', prob: [0.4, 0.05]},    
	    {obs: 'C', prob: [0.1, 0.4]},    
	    {obs: 'G', prob: [0.1, 0.4]},    
	    {obs: 'T', prob: [0.4, 0.05]}    
    ];

The initial probabilities of the hidden states will be a simple array.
Each element of the array has the same index of the corresponding hidden state in the hidden states array.

`let hiddenInit = [0.65, 0.35];`

In this example, 0.65 is the AT-rich probability, and the final 0.35 is the CG-rich probability.


To instantiate the Hidden Markov Model we pass the states, the observables and the initial probabilities as parameters of the HMM function.

`let HMModel = markov.HMM(hiddenStates, observables, hiddenInit);`

### Hidden Markov Model Bayes Theorem

To calculate the probability of a specific hidden state given an observable we can call the bayesTheorem() function and pass two parameters: the observable and the hidden state of which we want to know the probability.  

    let observation = 'A';
    let hiddenState = 'AT-rich';
    let bayesResult = HMModel.bayesTheorem(observation, hiddenState); //0.9369369369369369
  
### Hidden Markov Model Viterbi Algorithm

To calculate the most likely sequence of hidden states given a specific sequence of observables we can call the viterbiAlgorithm() function and pass it the observable sequence.

    let obSequence = ['A','T','C', 'G', 'C', 'G','T','C','A','T','C', 'G','T','C', 'G','T','C','C', 'G']; 
    let viterbiResult = HMModel.viterbiAlgorithm(obSequence);

The viterbiAlgorithm() returns an object with the following properties:

* states : the array of the hidden state sequence found

* prob : the resulting highest probability of the states sequence

* statesTrellis : an array of the trellis values of each state of the sequence. 

So, 

    console.log(viterbiResult.states) //[ 'AT-rich', 'AT-rich', 'AT-rich', 'AT-rich', 'CG-rich', 'CG-rich', ... ] 

### Hidden Markov Model properties

|Property | Description|
|------------ | -------------|
|states | Array of the names of the states|
|transMatrix | Array of arrays representing the transition probabilities|
|initialProb | Array of initial probabilities|
|observables | Array of the names of the observables|
|emissionMatrix | Array of arrays representing the emission probabilities|
|sequence | Sequence array of the input|
|sequenceArr | Index representation of the sequence array compared to the states array|


Example:

    console.log(HMModel.transMatrix) // [ [ 0.95, 0.05 ], [ 0.1, 0.9 ] ]