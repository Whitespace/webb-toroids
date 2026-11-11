"""Tell the user about toroids in a run that no other run has. Called by the queue worker.

Usage: queue_notify.py runs/NAME [--done]
A result is "new" if no other run under runs/ has one with the same group, genus, vertex, edge
and face counts, circumradius and pieces. Each new one is reported once (queue/state/notified.txt).
One notification per call, showing the simplest new toroid as a thumbnail.
On macOS it uses terminal-notifier if installed (thumbnail, click opens the picture), otherwise a
plain osascript notification. Elsewhere it prints."""
import csv, glob, os, platform, shutil, subprocess, sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
state = os.path.join(root, "queue", "state")


def rows(run):
    try:
        return list(csv.DictReader(open(os.path.join(run, "results.tsv")), delimiter="\t"))
    except OSError:
        return []


def key(r):
    return "|".join(r[k] for k in ("group", "genus", "V", "E", "F", "rcirc", "content"))


def notify(title, message, image=None):
    print("NOTIFY: %s: %s%s" % (title, message, " [%s]" % image if image else ""))
    if os.environ.get("NOTIFY", "1") != "1" or platform.system() != "Darwin":
        return
    tn = shutil.which("terminal-notifier")
    if tn:
        cmd = [tn, "-title", title, "-message", message, "-group", "webb-toroids", "-sound", "Glass"]
        if image:
            cmd += ["-contentImage", image, "-open", "file://" + image]
        subprocess.run(cmd, check=False)
    else:
        esc = lambda s: s.replace("\\", "\\\\").replace('"', '\\"')
        subprocess.run(["osascript", "-e", 'display notification "%s" with title "%s"' % (esc(message), esc(title))], check=False)


def short(content):
    for a, b in (("pentagonal-", "5-"), ("triangular-", "3-"), ("square-", "4-"), ("hexagonal-", "6-"), ("octagonal-", "8-"), ("decagonal-", "10-")):
        content = content.replace(a, b)
    return content.replace("+", " + ")


def main():
    run = os.path.normpath(sys.argv[1])
    done = "--done" in sys.argv
    os.makedirs(os.path.join(state, "thumbs"), exist_ok=True)
    npath = os.path.join(state, "notified.txt")
    notified = set(open(npath).read().split("\n")) if os.path.exists(npath) else set()
    known = set()
    for other in glob.glob(os.path.join(root, "runs", "*")):
        if os.path.abspath(other) != os.path.abspath(run):
            known |= {key(r) for r in rows(other)}
    mine = [r for r in rows(run) if r["class"] == "RAQT"]
    novel = [r for r in mine if key(r) not in known]
    fresh = [r for r in novel if key(r) not in notified]
    name = os.path.basename(run)
    if fresh:
        best = min(fresh, key=lambda r: int(r["F"]))
        png = os.path.join(state, "thumbs", "%s-%s.png" % (name, best["id"]))
        try:
            sys.path.insert(0, here)
            import render
            render.render(os.path.join(run, "off", "%06d.off" % int(best["id"])), png)
        except Exception as e:  # a failed thumbnail must not lose the notification
            print("thumbnail failed: %s" % e, file=sys.stderr)
            png = None
        title = "New toroid" if len(fresh) == 1 else "%d new toroids" % len(fresh)
        notify(title, "%s: genus %s, %s pieces, %s faces. %s" % (name, best["genus"], best["pieces"], best["F"], short(best["content"])), png)
        with open(npath, "a") as f:
            for r in fresh:
                f.write(key(r) + "\n")
    if done:
        notify("Search job finished", "%s: %d results, %d not in any other run" % (name, len(mine), len(novel)))


if __name__ == "__main__":
    main()
