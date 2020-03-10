import setuptools

if __name__ == '__main__':
    setuptools.setup(
        name='de_analysis_clean',
        version='0.1.0-dev',
        author='Erika Dudkin',
        author_email='erikadudkin@gmx.de',
        description='DE-Analysis for multi-patient groups',
        license='MIT License',
        url='url = https://github.com/erikadudki/de_analysis_clean',

        # Look for python code inside /src/
        packages=setuptools.find_packages(where='src'),

        # Assign the package-level folder ('') to be the one it finds in 'src'
        package_dir={'': 'src'},
        python_requires='>=3.6',

        install_requires=[
            'click'
        ]
    )