import java.util.concurrent.CountDownLatch;

public class SolveILPTask implements Runnable{
    private final Instance i;
    private final Tree T;
    private final CountDownLatch latch;
    private final ILPResultCollector coll;

    public ILPResult result = null;
    public ILPResult noHResult = null;
    public SolveILPTask(Instance i, Tree T, CountDownLatch latch, ILPResultCollector coll){
        this.i = i;
        this.T = T;
        this.latch = latch;
        this.coll = coll;
    }

    @Override
    public void run() {
        result = Calder.inferU(i, T);
        if (result != null){
            coll.addResult(result);
        } else {
            coll.addBadResult(T);
            noHResult = Calder.inferUwithoutH(i, T);
            if (noHResult != null)
                coll.addNoHResult(noHResult);
        }

        latch.countDown();
    }
}
