#include "Integrator.h"

#include "boost_for_export.h"
#include "globalDefs.h"
#include "cutils_func.h"
#include "DataSetUser.h"
#include "Fix.h"
#include "GPUArray.h"
#include "PythonOperation.h"
#include "WriteConfig.h"

using namespace std;




__global__ void zeroVectorPreserveW(float4 *xs, int n) {
    int idx = GETIDX();
    if (idx < n) {
        float w = xs[idx].w;
        xs[idx] = make_float4(0, 0, 0, w);
    }
}


void Integrator::stepInit(bool computeVirials)
{
    /*
    if (computeVirials) {
        //reset virials each turn
        state->gpd.virials.d_data.memset(0);
    }
    */

    for (Fix *f : state->fixes) {
        if (state->turn % f->applyEvery == 0) {
            f->stepInit();
        }
    }
    if (computeVirials) {
        //reset virials each turn
        state->gpd.virials.d_data.memset(0);
    }
}

void Integrator::stepFinal()
{
    for (Fix *f : state->fixes) {
        if (state->turn % f->applyEvery == 0) {
            f->stepFinal();
        }
    }
}








void Integrator::asyncOperations() {
    int turn = state->turn;

    // well, if I try to use a local state pointer, this segfaults. Need to
    // capture this instead.  Little confused
    auto writeAndPy = [this] (int64_t ts) {
        // have to set device in each thread
        state->devManager.setDevice(state->devManager.currentDevice, false);
        for (SHARED(WriteConfig) wc : state->writeConfigs) {
            if (not (ts % wc->writeEvery)) {
                wc->write(ts);
            }
        }
        for (SHARED(PythonOperation) po : state->pythonOperations) {
            if (not (ts % po->operateEvery)) {
                po->operate(ts);
            }
        }
    };
    bool needAsync = false;
    for (SHARED(WriteConfig) wc : state->writeConfigs) {
        if (not (turn % wc->writeEvery)) {
            needAsync = true;
            break;
        }
    }
    if (not needAsync) {
        for (SHARED(PythonOperation) po : state->pythonOperations) {
            if (not (turn % po->operateEvery)) {
                needAsync = true;
                break;
            }
        }
    }
    if (needAsync) {
        state->asyncHostOperation(writeAndPy);
    }
}



void Integrator::basicPreRunChecks() {
    if (state->devManager.prop.major < 3) {
        cout << "Device compute capability must be >= 3.0. Quitting" << endl;
        assert(state->devManager.prop.major >= 3);
    }
    if (state->rCut == RCUT_INIT) {
        cout << "rcut is not set" << endl;
        assert(state->rCut != RCUT_INIT);
    }
    if (state->is2d and state->periodic[2]) {
        cout << "2d system cannot be periodic is z dimension" << endl;
        assert(not (state->is2d and state->periodic[2]));
    }
    mdAssert(state->bounds.isInitialized(), "Bounds must be initialized");
    /*
    if (not state->bounds.isInitialized()) {
        cout << "Bounds not initialized" << endl;
        assert(state->bounds.isInitialized());
    }
    */

}


std::vector<bool> Integrator::basicPrepare(int numTurns) {
    std::cout << "Running for " << numTurns << " turns with timestep of " << state->dt << std::endl;
    int nAtoms = state->atoms.size();
    state->runningFor = numTurns;
    state->runInit = state->turn;
    state->nlistBuildCount = 0;
    state->prepareForRun();
    state->atomParams.guessAtomicNumbers();
    setActiveData();
    for (GPUArray *dat : activeData) {
        dat->dataToDevice();
    }
    std::vector<bool> prepared;
    for (Fix *f : state->fixes) {
        f->takeStateNThreadPerBlock(state->nThreadPerBlock);//grid will also have this value
        f->takeStateNThreadPerAtom(state->nThreadPerAtom);//grid will also have this value
        f->updateGroupTag();
        prepared.push_back(f->prepareForRun());
        f->setVirialTurnPrepare();
    }
    state->handleChargeOffloading();
    for (Fix *f : state->fixes) {
        f->setEvalWrapper(); //have to do this after prepare b/c pair calcs need evaluators from charge that have been updated with correct alpha or other coefficiants, and change calcs need to know that handoffs happened
    }
    state->gridGPU.periodicBoundaryConditions(-1, true);
    for (boost::shared_ptr<MD_ENGINE::DataSetUser> ds : state->dataManager.dataSets) {
        ds->prepareForRun(); //will also prepare those data sets' computers
        if (ds->requiresVirials()) {
            state->dataManager.addVirialTurn(ds->nextCompute, ds->requiresPerAtomVirials());
        }
    }

    return prepared;
}


void Integrator::basicFinish() {
    for (Fix *f : state->fixes) {
        f->postRun();
        f->hasAcceptedChargePairCalc = false;
        f->hasOffloadedChargePairCalc = false;
    }
    if (state->asyncData && state->asyncData->joinable()) {
        state->asyncData->join();
    }
    for (GPUArray *dat : activeData) {
        dat->dataToHost();
    }
    cudaDeviceSynchronize();
    state->downloadFromRun();
    state->finish();
}


void Integrator::setActiveData() {
    activeData = vector<GPUArray *>();
    activeData.push_back((GPUArray *) &state->gpd.ids);
    activeData.push_back((GPUArray *) &state->gpd.xs);
    activeData.push_back((GPUArray *) &state->gpd.vs);
    activeData.push_back((GPUArray *) &state->gpd.fs);
    activeData.push_back((GPUArray *) &state->gpd.idToIdxs);
    if (state->requiresCharges) {
        activeData.push_back((GPUArray *) &state->gpd.qs);
    }

    activeData.push_back((GPUArray *) &state->gpd.virials);
}


Integrator::Integrator(State *state_) : IntegratorUtil(state_) {
}


void Integrator::writeOutput() {
    for (SHARED(WriteConfig) wc : state->writeConfigs) {
        wc->write(state->turn);
    }
}


    /*
double Integrator::singlePointEngPythonAvg(string groupHandle) {
    GPUArrayGlobal<float> eng(2);
    eng.d_data.memset(0);
    basicPreRunChecks();
    basicPrepare(0);
    cudaDeviceSynchronize();

    singlePointEng();
    cudaDeviceSynchronize();
    uint32_t groupTag = state->groupTagFromHandle(groupHandle);
    int warpSize = state->devManager.prop.warpSize;
    accumulate_gpu_if<float, float, SumSingleIf, N_DATA_PER_THREAD> <<<NBLOCK(state->atoms.size() / (double) N_DATA_PER_THREAD), PERBLOCK, N_DATA_PER_THREAD*sizeof(float)*PERBLOCK>>>
        (
         eng.getDevData(), 
         state->gpd.perParticleEng.getDevData(),
         state->atoms.size(),
         warpSize,
         SumSingleIf(state->gpd.fs.getDevData(), groupTag)
        );
    eng.dataToHost();
    cudaDeviceSynchronize();
    CUT_CHECK_ERROR("Calculation of single point average energy failed");
    basicFinish();
    return eng.h_data[0] / *((int *)eng.h_data.data()+1);
}

boost::python::list Integrator::singlePointEngPythonPerParticle() {
    basicPrepare(0);
    singlePointEng();
    state->gpd.perParticleEng.dataToHost();
    state->gpd.ids.dataToHost();
    cudaDeviceSynchronize();
    CUT_CHECK_ERROR("Calculation of single point per-particle energy failed");
    vector<float> &engs = state->gpd.perParticleEng.h_data;
    vector<uint> &ids = state->gpd.ids.h_data;
    vector<int> &idToIdxsOnCopy = state->gpd.idToIdxsOnCopy;
    vector<double> sortedEngs(ids.size());

    for (int i=0, ii=state->atoms.size(); i<ii; i++) {
        int id = ids[i];
        int idxWriteTo = idToIdxsOnCopy[id];
        sortedEngs[idxWriteTo] = eng:tas[i];
    }
    boost::python::list asPy(sortedEngs);
    basicFinish();
    return asPy;

}

    */


double Integrator::tune() {
    return 0;
    
    auto startTune = std::chrono::high_resolution_clock::now();

    std::vector<int> threadPerBlocks= {32, 64, 128, 256};
    std::vector<int> threadPerAtoms = {1, 2, 4, 8, 16, 32};
    //if we have run for a while, guess the fraction of turns where we build nlists, else guess that we build every other time
    double nlistBuildFrac = state->turn > state->runInit ? (double) state->nlistBuildCount/(state->turn-state->runInit) : 1/(2.0 * state->periodicInterval);
    int nForceEvals = 30;
    //estimating how many times we should build nlists for good estimate
    int nNlistBuilds = round(nForceEvals * nlistBuildFrac);

    auto setParams = [this] (int ntpb, int ntpa) {
        state->nThreadPerBlock = ntpb;
        state->nThreadPerAtom = ntpa;
        state->gridGPU.nThreadPerBlock(ntpb);
        state->gridGPU.nThreadPerAtom(ntpa);
        state->gridGPU.initArraysTune();
        for (auto fix : state->fixes) {
            fix->takeStateNThreadPerBlock(ntpb);
            fix->takeStateNThreadPerAtom(ntpa);
        }
    };

    vector<vector<double> > times;
    //REMEMBER TO MAKE COPY OF FORCES AND SET THEM BACK AFTER THIS;
    for (int i=0; i<threadPerBlocks.size(); i++) {
        vector<double> timesWithBlock;
        for (int j=0; j<threadPerAtoms.size(); j++) {
            int threadPerBlock = threadPerBlocks[i];
            int threadPerAtom = threadPerAtoms[j];
            setParams(threadPerBlock, threadPerAtom);
            state->gridGPU.periodicBoundaryConditions(-1, true);
            state->nlistBuildCount--;
            cudaDeviceSynchronize();
            auto start = std::chrono::high_resolution_clock::now();
            for (int k=0; k<nNlistBuilds; k++) {
                state->gridGPU.periodicBoundaryConditions(-1, true);
                state->nlistBuildCount--;
            }
            for (int k=0; k<nForceEvals; k++) {
                force(false);
            }
            cudaDeviceSynchronize();
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start;
            timesWithBlock.push_back(duration.count());
        }
        times.push_back(timesWithBlock);
    }

    double bestTime = times[0][0];
    int bestNTPB, bestNTPA;
    bestNTPB = 0;
    bestNTPA = 0;
    
    
    
    for (int i=0; i<threadPerBlocks.size(); i++) {
        for (int j=0; j<threadPerAtoms.size(); j++) {
            if (times[i][j] < bestTime) {
                bestNTPB=threadPerBlocks[i];
                bestNTPA=threadPerAtoms[j];
                bestTime = times[i][j];
            }
        }
    }
    setParams(bestNTPB, bestNTPA);
    state->gridGPU.periodicBoundaryConditions(-1, true);
    state->nlistBuildCount--;
    //then pick best one


    //zero forces that you calculated here
    zeroVectorPreserveW<<<NBLOCKVAR(state->atoms.size(), state->nThreadPerBlock), state->nThreadPerBlock>>>(state->gpd.fs.getDevData(), state->atoms.size());
    auto endTune = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> durationTune = endTune - startTune;
    return durationTune.count();
}

void export_Integrator() {
    boost::python::class_<Integrator,
                          boost::noncopyable> (
        "Integrator"
    )
    .def("writeOutput", &Integrator::writeOutput)
   // .def("energyAverage", &Integrator::singlePointEngPythonAvg,
    //        (boost::python::arg("groupHandle")="all")
    //    )
    //.def("energyPerParticle", &Integrator::singlePointEngPythonPerParticle);
    //.def("run", &Integrator::run)
    ;
}

