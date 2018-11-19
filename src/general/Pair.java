package general;

public class Pair<E> {
	public final E a;
	public final E b;

	public Pair(E a, E b) {
		this.a = a;
		this.b = b;
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof Pair)) {
			return false;
		}

		@SuppressWarnings("unchecked")
		Pair<E> other = (Pair<E>) o;
		if (this.a.equals(other.a) && this.b.equals(other.b)) {
			return true;
		} else {
			return false;
		}
	}

	@Override
	public int hashCode() {
		return a.hashCode() + (31 * b.hashCode());
	}
}