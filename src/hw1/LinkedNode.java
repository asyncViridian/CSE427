package hw1;

import general.Pair;

public class LinkedNode {
    public final Pair<Integer> coords;
    public final int value;
    public final LinkedNode prev;

    public LinkedNode(Pair<Integer> coords, int value, LinkedNode prev) {
        this.coords = coords;
        this.value = value;
        this.prev = prev;
    }
}