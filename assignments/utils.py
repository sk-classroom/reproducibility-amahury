from matplotlib.colors import ListedColormap
import numpy as np
import matplotlib.pyplot as plt


def plot_decision_regions(X, y, classifier, test_idx=None, resolution=1000):
    """
    Plot decision regions of a classifier.

    Parameters:
    X : array-like, shape = [n_samples, n_features]
        Feature Matrix.
    y : array-like, shape = [n_samples]
        True class labels.
    classifier : Classifier object. Must have a .predict method.
    test_idx : array-like, shape = [n_samples]
        Indices of test set examples.
    resolution : float, optional (default: 0.02)
        The resolution of the grid for plotting decision regions.

    Returns:
    None
    """
    # setup marker generator and color map
    markers = ("o", "s", "^", "v", "<")
    colors = ("red", "blue", "lightgreen", "gray", "cyan")
    cmap = ListedColormap(colors[: len(np.unique(y))])

    # plot the decision surface
    x1_min, x1_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    x2_min, x2_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    resolution_x1 = (x1_max - x1_min) / resolution
    resolution_x2 = (x2_max - x2_min) / resolution
    xx1, xx2 = np.meshgrid(
        np.arange(x1_min, x1_max, resolution_x1),
        np.arange(x2_min, x2_max, resolution_x2),
    )
    lab = classifier.predict(np.array([xx1.ravel(), xx2.ravel()]).T)
    lab = np.unique(lab, return_inverse=True)[1]
    lab = lab.reshape(xx1.shape)
    plt.contourf(xx1, xx2, lab, alpha=0.3, cmap=cmap)
    plt.xlim(xx1.min(), xx1.max())
    plt.ylim(xx2.min(), xx2.max())

    # plot class examples
    for idx, cl in enumerate(np.unique(y)):
        plt.scatter(
            x=X[y == cl, 0],
            y=X[y == cl, 1],
            alpha=0.8,
            c=colors[idx],
            marker=markers[idx],
            label=f"Class {cl}",
            edgecolor="black",
        )

    # highlight test examples
    if test_idx:
        # plot all examples
        X_test, y_test = X[test_idx, :], y[test_idx]

        plt.scatter(
            X_test[:, 0],
            X_test[:, 1],
            c="none",
            edgecolor="black",
            alpha=1.0,
            linewidth=1,
            marker="o",
            s=100,
            label="Test set",
        )
    plt.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        shadow=True,
        ncol=2,
        frameon=False,
    )
