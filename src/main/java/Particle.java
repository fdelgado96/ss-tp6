import java.util.Locale;

public class Particle {

    public double x, y, vx, vy, fx, fy, m = 80, prevX, prevY, r, fn, perimeter;
    public boolean isOut = false;

    private boolean initialized = false;
    private int id;

    public Particle(int id, double x, double y, double r, boolean thicc_mode) {
        this.x = x;
        this.y = y;
        this.r = r;
        if  (thicc_mode) {
            this.m = r * 700 - 105;
        }
        this.id = id;
        this.perimeter = Math.PI * 2 * r;
    }


    public double getOverlap(Particle other) {
        double overlap = r + other.r - centerDistance(other);
        return overlap >= 0 ? overlap : 0;
    }

    public double centerDistance(Particle other){
        return Math.sqrt(Math.pow(other.x - x, 2) + Math.pow(other.y - y, 2));
    }

    public void clearForces() {
        fn = 0;
        fx = 0;
        fy = 0;
    }

    public void clearVelocities() {
        vx = 0;
        vy = 0;
        prevX = x;
        prevY = y;
        initialized = false;
    }

    public double enx(Particle other) {
        return (other.x - x)/centerDistance(other);
    }

    public double centerDistanceToPoint(double x, double y){
        return Math.sqrt(Math.pow(x - this.x, 2) + Math.pow(y - this.y, 2));
    }

    public double eTargetX(double x, double y) {
        return (x - this.x)/centerDistanceToPoint(x, y);
    }
    public double eTargetY(double x, double y) {
        return (y - this.y)/centerDistanceToPoint(x, y);
    }

    public double eny(Particle other) {
        return (other.y - y)/centerDistance(other);
    }

    public double kineticEnergy() {
        return 0.5 * m * (vx*vx + vy*vy);
    }

    public void move(double dt) {
        if (!initialized) {
            prevX = x - vx * dt;
            prevY = y - vy * dt;
            initialized = true;
        }

        double dt2  = dt*dt;
        double nextX = 2*x - prevX + fx * dt2 / m;
        double nextY = 2*y - prevY + fy * dt2 / m;

        prevX = x;
        prevY = y;
        x = nextX;
        y = nextY;
        vx = (x - prevX) / dt;
        vy = (y - prevY) / dt;
    }


    public double getNormalRelVel(Particle p2) {
        return ((vx - p2.vx) * enx(p2)) + ((vy - p2.vy) * eny(p2));
    }

    @Override
    public String toString() {
        return String.format(Locale.ENGLISH, "%d %.12f %.12f %.4f %.8f", id, x, y, r, fn/perimeter);
    }


    public double getTangentRelVel(Particle p2) {
        return ((vx - p2.vx) * -eny(p2)) + ((vy - p2.vy) * enx(p2));
    }
}
