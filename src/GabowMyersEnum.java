public class GabowMyersEnum extends GabowMyers {

    public GabowMyersEnum(Graph G, int nSamples){
        this.myG = G;
        this.nSamples = nSamples;
    }

    @Override
    protected void output(Tree T, int root) {
        L = new Tree(T, root, nSamples);
        if(L.isValid()){
            results.add(L);
        }
    }
}
